# =============================================================================
# NB_Comparisons.R
# Mixtures of shifted negative binomial distributions: one sampler, five updates of the
# component parameters, for the comparison asked for in Referee 1, comment 4b.
#
# Everything except the update of the negative binomial parameters of each non-empty component
# is shared: the allocation step (lambda integrated out), the joint update of the shift and the
# mean, the drawing of empty components from their prior, the sparse prior on the weights
# (alpha by random walk on log alpha), the weights, the adaptation rule for random walks and
# the stored output. Each algorithm is selected with `algorithm`; with the same data, starting
# values and seed, the runs differ only in that update.
#
# Model: y_i - kappa_j ~ NB(mu_j, r_j) given z_i = j, with mean mu_j and variance mu_j + mu_j^2 / r_j.
# Our prior: mu_j ~ Unif(0, M], rtilde_j = 1 / (1 + r_j) ~ Unif(0, 1/2), kappa_j ~ DiscUnif{0, ..., M},
# pi ~ Dirichlet(alpha), alpha ~ Gamma(1, rate 200), with M = max(y).
#
# algorithm =
#  "ours"    : rtilde by random walk with lambda integrated out, then lambda and mu given lambda
#              (mu is drawn exactly from a truncated Inverse-Gamma). Our prior.
#  "mur_rw"  : rtilde by random walk as in "ours", then log(mu) by random walk. Our prior.
#  "pr_rw"   : log(r) by random walk given p, then logit(p) by random walk given r, with
#              p = r / (r + mu). Our prior, transformed to (r, p).
#  "zhou"    : the Gibbs update of Zhou et al. (2012) and Zhou and Carin (2012), NegBino_Gibbs.m,
#              lines 79-81, under Zhou's priors q ~ Beta(a0, b0), r ~ Gamma(e0, rate f0).
#  "zhou_mh" : the same three draws, used as Metropolis-Hastings proposals with our prior, so that
#              it targets the same posterior as "ours", "mur_rw" and "pr_rw".
#  "sp"      : the mixture of shifted Poisson distributions (Cross et al., 2024), i.e. our model with
#              r_j = Inf (stored as rtilde_j = 0). Same prior on mu_j and kappa_j, same allocation, joint
#              kappa step and empty components; mu_j given the allocation is drawn exactly from
#              Gamma(sum(y - kappa) + 1, rate n_j) truncated to (0, M].
# Zhou's probability q is our 1 - p: his pmf is proportional to q^x (1 - q)^r, with mean r q / (1 - q).
#
# Random walks are adapted in batches of 50 iterations during burn-in only (target acceptance
# 0.3 to 0.4) and fixed afterwards. Run time is CPU time of the sampling loop only.
# Base R only.
# =============================================================================

ZHOU_HP <- c(a0 = 1, b0 = 1, e0 = 1, f0 = 0.1)   # Zhou's priors as used here: q uniform, r exponential, mean 10

nb_ll <- function(ya, na, r, mu) sum(na * dnbinom(ya, size = r, mu = mu, log = TRUE))

# ---- Shared: joint update of kappa_j and mu_j (Section 3.1) -----------------------------------
# kappa_j is drawn given m_j = kappa_j + mu_j and r_j, then mu_j = m_j - kappa_j. The change of
# coordinates has a unit Jacobian. Under our prior mu is uniform, so only the likelihood enters;
# under Zhou's prior the density of mu given r implied by q ~ Beta(a0, b0), q = mu / (mu + r), enters,
# r / (mu + r)^2 * Beta(q), and mu has no upper bound.
draw_kappa <- function(ua, na, kappa, mu, r, M, zhou_prior = FALSE, hp = ZHOU_HP) {
  m <- kappa + mu
  grid <- 0:min(min(ua), M)
  grid <- if (zhou_prior) grid[(m - grid) > 0] else grid[(m - grid) > 0 & (m - grid) <= M]
  lp <- vapply(grid, function(k) {
    mk <- m - k
    out <- sum(na * dnbinom(ua - k, size = r, mu = mk, log = TRUE))
    if (zhou_prior) out <- out + dbeta(mk / (mk + r), hp["a0"], hp["b0"], log = TRUE) + log(r) - 2 * log(mk + r)
    out
  }, numeric(1))
  k_new <- if (length(grid) == 1) grid else grid[sample.int(length(grid), 1, prob = exp(lp - max(lp)))]
  list(kappa = k_new, mu = m - k_new)
}

# ---- "ours" and "mur_rw" ------------------------------------------------------------------
draw_rtilde <- function(ya, na, mu, rtilde, sd_rw) {
  cand <- rnorm(1, rtilde, sd_rw)
  if (cand <= 0 || cand >= 0.5) return(list(value = rtilde, acc = 0))
  if (log(runif(1)) < nb_ll(ya, na, (1 - cand) / cand, mu) - nb_ll(ya, na, (1 - rtilde) / rtilde, mu))
    list(value = cand, acc = 1) else list(value = rtilde, acc = 0)
}
# lambda, then mu given lambda: 1/mu is Gamma(r n - 1, rate r sum(lambda)) truncated to [1/M, Inf).
# Only the sum of the lambda_i is needed; it is one Gamma draw.
draw_mu <- function(ya, na, mu, r, M) {
  n <- sum(na)
  sum_lambda <- rgamma(1, shape = sum(na * ya) + n * r, scale = mu / (mu + r))
  a <- n * r - 1; b <- r * sum_lambda
  upper <- pgamma(1 / M, a, rate = b, lower.tail = FALSE)
  if (upper < 1e-12) return(M)
  1 / qgamma(runif(1) * upper, a, rate = b, lower.tail = FALSE)
}
# Shifted Poisson: with mu uniform on (0, M], mu given the data is Gamma(sum(y) + 1, rate n) truncated to (0, M]
draw_mu_pois <- function(ya, na, M) {
  a <- sum(na * ya) + 1; b <- sum(na)
  upper <- pgamma(M, a, rate = b)
  if (upper < 1e-300) return(M)
  qgamma(runif(1) * upper, a, rate = b)
}
draw_logmu <- function(ya, na, mu, r, M, sd_rw) {          # uniform prior on mu; Jacobian mu
  cand <- mu * exp(sd_rw * rnorm(1))
  if (cand > M) return(list(value = mu, acc = 0))
  if (log(runif(1)) < nb_ll(ya, na, r, cand) - nb_ll(ya, na, r, mu) + log(cand) - log(mu))
    list(value = cand, acc = 1) else list(value = mu, acc = 0)
}

# ---- "pr_rw": our prior in (r, p), p = r / (r + mu), density proportional to (1 + r)^-2 r / p^2 ----
lprior_rp <- function(r, p) -2 * log(1 + r) + log(r) - 2 * log(p)
draw_logr_given_p <- function(ya, na, r, p, M, sd_rw) {    # Jacobian r
  cand <- r * exp(sd_rw * rnorm(1)); mu_c <- cand * (1 - p) / p
  if (cand <= 1 || mu_c > M) return(list(value = r, acc = 0))
  la <- nb_ll(ya, na, cand, mu_c) - nb_ll(ya, na, r, r * (1 - p) / p) + lprior_rp(cand, p) - lprior_rp(r, p) + log(cand) - log(r)
  if (log(runif(1)) < la) list(value = cand, acc = 1) else list(value = r, acc = 0)
}
draw_logitp_given_r <- function(ya, na, r, p, M, sd_rw) {  # Jacobian p (1 - p)
  cand <- plogis(qlogis(p) + sd_rw * rnorm(1)); mu_c <- r * (1 - cand) / cand
  if (mu_c > M || mu_c <= 0) return(list(value = p, acc = 0))
  la <- nb_ll(ya, na, r, mu_c) - nb_ll(ya, na, r, r * (1 - p) / p) + lprior_rp(r, cand) - lprior_rp(r, p) +
    log(cand * (1 - cand)) - log(p * (1 - p))
  if (log(runif(1)) < la) list(value = cand, acc = 1) else list(value = p, acc = 0)
}

# ---- "zhou" and "zhou_mh" -------------------------------------------------------------------
# Zhou's CRT_sum draws, for each observation x_i, sum_{k = 1}^{x_i} Bernoulli(r / (r + k - 1)).
# With unique values ya and counts na this is sum_k Binomial(number of observations with ya >= k,
# r / (r + k - 1)), which has the same distribution. The k = 1 probability is 1; pmin removes rounding.
crt_sum <- function(ya, na, r) {
  K <- max(ya); if (K == 0) return(0)
  sum(rbinom(K, vapply(1:K, function(k) sum(na[ya >= k]), numeric(1)), pmin(1, r / (r + (1:K) - 1))))
}
# Our prior in Zhou's coordinates (r, q), q = 1 - p, mu = r q / (1 - q): proportional to
# (1 + r)^-2 r / (1 - q)^2 on r > 1 and mu <= M. Zhou's prior: Beta(q; a0, b0) Gamma(r; e0, f0).
lprior_ours_rq <- function(r, q, M) {
  if (r <= 1 || r * q / (1 - q) > M) return(-Inf)
  -2 * log(1 + r) + log(r) - 2 * log(1 - q)
}
lprior_zhou_rq <- function(r, q, hp) dbeta(q, hp["a0"], hp["b0"], log = TRUE) + dgamma(r, hp["e0"], rate = hp["f0"], log = TRUE)
zhou_update <- function(ya, na, r, q, M, hp, mh) {
  N <- sum(na); acc <- c(1, 1)
  # 1. q from Beta(a0 + sum x, b0 + N r): exact under Zhou's prior; independence proposal under ours
  q_c <- rbeta(1, hp["a0"] + sum(na * ya), hp["b0"] + N * r)
  if (mh) {
    lw <- function(qq) lprior_ours_rq(r, qq, M) - dbeta(qq, hp["a0"], hp["b0"], log = TRUE)
    if (log(runif(1)) < lw(q_c) - lw(q)) q <- q_c else acc[1] <- 0
  } else q <- q_c
  # 2. L ~ CRT given r (exact under either prior); 3. r from Gamma(e0 + L, rate f0 - N log(1 - q))
  L <- crt_sum(ya, na, r)
  r_c <- rgamma(1, hp["e0"] + L, rate = hp["f0"] - N * log(max(1 - q, .Machine$double.xmin)))
  if (mh) {
    lw <- function(rr) lprior_ours_rq(rr, q, M) - dgamma(rr, hp["e0"], rate = hp["f0"], log = TRUE)
    if (log(runif(1)) < lw(r_c) - lw(r)) r <- r_c else acc[2] <- 0
  } else r <- r_c
  list(r = unname(r), q = unname(q), acc = acc)
}

# ---- Shared: alpha, by random walk on log(alpha) with the weights integrated out -------------
draw_alpha <- function(alpha, nj, step = 1, shape0 = 1, rate0 = 200) {
  J <- length(nj); n <- sum(nj)
  lt <- function(a) (shape0 - 1) * log(a) - rate0 * a + lgamma(J * a) - lgamma(n + J * a) + sum(lgamma(nj + a) - lgamma(a))
  cand <- exp(log(alpha) + step * rnorm(1))
  if (log(runif(1)) < lt(cand) - lt(alpha) + log(cand) - log(alpha)) cand else alpha
}

# ---- The sampler ----------------------------------------------------------------------------
# init: list(kappa, rtilde) of length J (optional). seed: set at the start, so equal seeds give
# equal starting allocations across algorithms. Returns the kept draws (every `thin`-th), acceptance
# rates of the two component updates, and CPU seconds of the sampling loop.
nb_sampler <- function(y, J = 10, burnin = 5000, nkeep = 10000, algorithm = "ours", thin = 1,
                       sparse = (J > 1), alpha = if (sparse) 0.005 else 1, shifted = TRUE,
                       init = NULL, hp = ZHOU_HP, seed = NULL) {
  stopifnot(algorithm %in% c("ours", "mur_rw", "pr_rw", "zhou", "zhou_mh", "sp"))
  sp <- algorithm == "sp"
  if (!is.null(seed)) set.seed(seed)
  zp <- algorithm == "zhou"                                  # Zhou's prior; all others use ours
  tab <- table(y); u <- as.numeric(names(tab)); cnt <- as.numeric(tab); m_u <- length(u); M <- max(y)
  kappa <- rep(0, J); rtilde <- rep(1 / 11, J); mu <- rep(max(mean(y), 1), J)
  if (!is.null(init)) { rtilde <- init$rtilde; if (shifted) kappa <- init$kappa }
  if (sp) rtilde <- rep(0, J)                                # r = Inf: Poisson components
  w <- rep(1 / J, J)
  sd1 <- rep(if (algorithm == "pr_rw") 0.3 else 0.1, J); sd2 <- rep(0.1, J)
  nuj <- matrix(cnt, m_u, 1)
  if (J > 1) nuj <- t(sapply(cnt, function(cc) rmultinom(1, cc, rep(1 / J, J))))
  nst <- floor(nkeep / thin)
  store <- list(mu = matrix(NA, nst, J), rtilde = matrix(NA, nst, J), kappa = matrix(NA, nst, J),
                w = matrix(NA, nst, J), nj = matrix(NA, nst, J), alpha = rep(NA, nst))
  b1 <- b2 <- n1 <- n2 <- k1 <- k2 <- kn1 <- kn2 <- rep(0, J); n_batch <- 0
  t0 <- proc.time()

  for (it in 1:(burnin + nkeep)) {
    # Allocation, lambda integrated out
    if (J > 1) {
      lp <- matrix(-Inf, m_u, J)
      for (j in 1:J) {
        id <- u >= kappa[j]
        lp[id, j] <- log(w[j]) + dnbinom(u[id] - kappa[j], size = (1 - rtilde[j]) / rtilde[j], mu = mu[j], log = TRUE)
      }
      pr <- exp(lp - apply(lp, 1, max))
      for (i in 1:m_u) nuj[i, ] <- rmultinom(1, cnt[i], pr[i, ])
    }
    nj <- colSums(nuj)

    for (j in 1:J) {
      if (nj[j] == 0) {                                      # empty component: draw from the prior
        if (shifted) kappa[j] <- sample.int(M + 1, 1) - 1
        if (zp) { q <- rbeta(1, hp["a0"], hp["b0"]); r <- rgamma(1, hp["e0"], rate = hp["f0"])
                  rtilde[j] <- 1 / (1 + r); mu[j] <- max(r * q / (1 - q), 1e-8) }
        else { mu[j] <- M * (1 - runif(1)); rtilde[j] <- if (sp) 0 else 0.5 * (1 - runif(1)) }
        next
      }
      act <- nuj[, j] > 0; ua <- u[act]; na <- nuj[act, j]
      r <- (1 - rtilde[j]) / rtilde[j]
      if (shifted) { s <- draw_kappa(ua, na, kappa[j], mu[j], r, M, zp, hp); kappa[j] <- s$kappa; mu[j] <- s$mu }
      ya <- ua - kappa[j]
      # ---- the only part that differs between algorithms ----
      if (algorithm == "ours") {
        s <- draw_rtilde(ya, na, mu[j], rtilde[j], sd1[j]); rtilde[j] <- s$value; a1 <- s$acc
        mu[j] <- draw_mu(ya, na, mu[j], (1 - rtilde[j]) / rtilde[j], M); a2 <- NA
      } else if (algorithm == "mur_rw") {
        s <- draw_rtilde(ya, na, mu[j], rtilde[j], sd1[j]); rtilde[j] <- s$value; a1 <- s$acc
        s <- draw_logmu(ya, na, mu[j], (1 - rtilde[j]) / rtilde[j], M, sd2[j]); mu[j] <- s$value; a2 <- s$acc
      } else if (algorithm == "pr_rw") {
        p <- r / (r + mu[j])
        s <- draw_logr_given_p(ya, na, r, p, M, sd1[j]); r <- s$value; a1 <- s$acc
        s <- draw_logitp_given_r(ya, na, r, p, M, sd2[j]); p <- s$value; a2 <- s$acc
        rtilde[j] <- 1 / (1 + r); mu[j] <- r * (1 - p) / p
      } else if (sp) {
        mu[j] <- draw_mu_pois(ya, na, M); a1 <- NA; a2 <- NA
      } else {                                               # "zhou" or "zhou_mh"
        s <- zhou_update(ya, na, r, mu[j] / (mu[j] + r), M, hp, mh = algorithm == "zhou_mh")
        rtilde[j] <- 1 / (1 + s$r); mu[j] <- max(s$r * s$q / (1 - s$q), 1e-8)
        a1 <- if (algorithm == "zhou_mh") s$acc[2] else NA; a2 <- if (algorithm == "zhou_mh") s$acc[1] else NA
      }
      # ---------------------------------------------------------
      if (!is.na(a1)) { b1[j] <- b1[j] + a1; n1[j] <- n1[j] + 1; if (it > burnin) { k1[j] <- k1[j] + a1; kn1[j] <- kn1[j] + 1 } }
      if (!is.na(a2)) { b2[j] <- b2[j] + a2; n2[j] <- n2[j] + 1; if (it > burnin) { k2[j] <- k2[j] + a2; kn2[j] <- kn2[j] + 1 } }
    }
    if (sparse) alpha <- draw_alpha(alpha, nj)
    g <- rgamma(J, nj + alpha, 1); if (!all(is.finite(g)) || sum(g) == 0) g <- nj + 1e-10
    w <- g / sum(g)

    # Adaptation of random walk scales, burn-in only (zhou and zhou_mh have none to adapt)
    if (it <= burnin && it %% 50 == 0) {
      n_batch <- n_batch + 1; dl <- min(0.1, 1 / sqrt(n_batch))
      if (algorithm %in% c("ours", "mur_rw", "pr_rw")) for (j in 1:J) {
        if (n1[j] > 0) { rt <- b1[j] / n1[j]; if (rt > 0.4) sd1[j] <- sd1[j] * exp(dl) else if (rt < 0.3) sd1[j] <- sd1[j] * exp(-dl) }
        if (n2[j] > 0) { rt <- b2[j] / n2[j]; if (rt > 0.4) sd2[j] <- sd2[j] * exp(dl) else if (rt < 0.3) sd2[j] <- sd2[j] * exp(-dl) }
      }
      b1[] <- 0; n1[] <- 0; b2[] <- 0; n2[] <- 0
    }
    if (it > burnin && (it - burnin) %% thin == 0) {
      i <- (it - burnin) / thin
      store$mu[i, ] <- mu; store$rtilde[i, ] <- rtilde; store$kappa[i, ] <- kappa
      store$w[i, ] <- w; store$nj[i, ] <- nj; store$alpha[i] <- alpha
    }
  }
  tt <- proc.time() - t0
  store$cpu <- unname(tt["user.self"] + tt["sys.self"])
  store$acc1 <- ifelse(kn1 > 0, k1 / pmax(kn1, 1), NA); store$acc2 <- ifelse(kn2 > 0, k2 / pmax(kn2, 1), NA)
  store$algorithm <- algorithm
  store$sd1 <- sd1; store$sd2 <- sd2                         # random-walk scales after the burn-in adaptation
  store
}

# ---- Summaries of a run (label-invariant) ---------------------------------------------------
# Number and location of modes for one draw (Section 3.2): non-empty components only.
nb_modes <- function(mu, rtilde, kappa, w, nj) {
  act <- nj > 0; mu <- mu[act]; rtilde <- rtilde[act]; kappa <- kappa[act]; w <- w[act] / sum(w[act])
  r <- (1 - rtilde) / rtilde
  top <- max(kappa + pmax(0, floor(ifelse(is.infinite(r), mu, (r - 1) * mu / r)))) + 1   # r = Inf: Poisson
  yg <- 0:top; pm <- rep(0, length(yg))
  for (j in seq_along(r)) { id <- yg >= kappa[j]; pm[id] <- pm[id] + w[j] * dnbinom(yg[id] - kappa[j], size = r[j], mu = mu[j]) }
  s <- sign(diff(c(0, pm))); nz <- which(s != 0); L <- 0
  for (t in seq_along(nz)) if (s[nz[t]] == 1 && (t == length(nz) || s[nz[t + 1]] == -1)) L <- L + 1
  L
}
# Log-likelihood of the mixture at each stored draw
nb_loglik <- function(o, y) {
  tab <- table(y); u <- as.numeric(names(tab)); cnt <- as.numeric(tab)
  vapply(seq_len(nrow(o$mu)), function(i) {
    act <- which(o$nj[i, ] > 0); w <- o$w[i, act] / sum(o$w[i, act]); pu <- rep(0, length(u))
    for (k in seq_along(act)) { j <- act[k]; id <- u >= o$kappa[i, j]
      pu[id] <- pu[id] + w[k] * dnbinom(u[id] - o$kappa[i, j], size = (1 - o$rtilde[i, j]) / o$rtilde[i, j], mu = o$mu[i, j]) }
    sum(cnt * log(pu))
  }, numeric(1))
}
# Effective sample size (Geyer's initial positive sequence) and inefficiency factor = n / ESS
nb_ess <- function(x) {
  n <- length(x); if (any(!is.finite(x)) || var(x) == 0) return(NA_real_)   # e.g. dispersion under "sp" (r = Inf)
  a <- acf(x, lag.max = min(n - 1, 1000), plot = FALSE)$acf[, 1, 1]
  m <- floor((length(a) - 1) / 2); g <- a[seq(1, 2 * m - 1, by = 2)] + a[seq(2, 2 * m, by = 2)]
  s <- which(g <= 0)[1]; if (!is.na(s)) g <- g[seq_len(s - 1)]
  if (!length(g)) return(n); min(n, n / (2 * sum(g) - 1))
}
nb_summary <- function(o, y) {
  ll <- nb_loglik(o, y)
  L <- vapply(seq_len(nrow(o$mu)), function(i) nb_modes(o$mu[i, ], o$rtilde[i, ], o$kappa[i, ], o$w[i, ], o$nj[i, ]), numeric(1))
  act <- o$nj > 0; ww <- o$w * act; ww <- ww / rowSums(ww)
  disp <- rowSums(ww * log((1 - o$rtilde) / o$rtilde))      # sum_j pi_j log r_j over non-empty components
  K <- length(ll); e <- c(loglik = nb_ess(ll), modes = nb_ess(L), dispersion = nb_ess(disp))
  list(algorithm = o$algorithm, cpu = o$cpu, loglik = ll, L = L, disp = disp,
       ess = e, if_ = K / e, ess_per_cpu_sec = e / o$cpu, pr_L = table(L) / K,
       acc = c(median(o$acc1, na.rm = TRUE), median(o$acc2, na.rm = TRUE)))
}

# ---- Run all algorithms on the same data, the same starting values and the same seed ---------
nb_compare <- function(y, algorithms = c("ours", "mur_rw", "pr_rw", "zhou", "zhou_mh"), seed = 1, ...) {
  lapply(setNames(algorithms, algorithms), function(a) nb_summary(nb_sampler(y, algorithm = a, seed = seed, ...), y))
}

# ---- Example (runs only with Rscript NB_Comparisons.R) ----------------------------------------
if (sys.nframe() == 0) {
  set.seed(4001)
  z <- sample.int(2, 2000, TRUE, prob = c(.6, .4))
  y <- c(0, 10)[z] + rnbinom(2000, size = c(5, 3)[z], mu = c(2, 5)[z])   # G4: shift of 10 in component 2
  res <- nb_compare(y, burnin = 2000, nkeep = 4000, thin = 5, seed = 7)
  for (a in names(res)) with(res[[a]], cat(sprintf("%-8s CPU %5.1f s | IF log-lik %5.1f modes %5.1f dispersion %6.1f | Pr(L = 2) %.2f\n",
    a, cpu, if_[1], if_[2], if_[3], if ("2" %in% names(pr_L)) pr_L[["2"]] else 0)))
}
