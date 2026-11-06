# -----------------------------------------------------------------------------
# Finite mixture of shifted negative binomial distributions: MCMC sampler and
# mode inference, following Sections 3.1 and 3.2 of the paper.
#
# Model:  y_i - kappa_j ~ NB(mu_j, rtilde_j) if z_ij = 1, with r_j = (1 - rtilde_j) / rtilde_j
# Priors: mu_j ~ Unif(0, M], rtilde_j ~ Unif(0, 1/2), kappa_j ~ DiscUnif(0, M),
#         pi ~ Dirichlet(alpha, ..., alpha), with M = max(y), and either alpha fixed or
#         alpha ~ Gamma(1, 200) (shape, rate).
#
# Blocks, in order: z, kappa, rtilde, lambda, mu, alpha, pi.
#
# joint_kappa chooses the update of kappa:
#   FALSE: kappa_j given mu_j, as Section 3.1 and Appendix C currently state.
#   TRUE:  kappa_j given m_j = kappa_j + mu_j, then mu_j = m_j - kappa_j. The change of
#          coordinates has a unit Jacobian and the priors are uniform, so this is a Gibbs
#          step in (kappa_j, m_j). The grid drops any kappa with mu_j outside (0, M].
#          In our checks, FALSE could leave kappa_j stuck, because given mu_j the data pin
#          kappa_j close to the sample mean minus mu_j. The argument needs checking by Lennart.
#
# Sum of lambda: only the sum of the lambda_i in a component is needed to draw mu_j. The lambda_i
# are independent Gamma with a common scale, so the sum is drawn directly as one Gamma variable.
#
# Base R only. Run with:  Rscript snb_mixture_sampler.R
# -----------------------------------------------------------------------------

# One draw of the sampler's kappa step ---------------------------------------
draw_kappa <- function(ua, na, kappa, mu, r, M, joint_kappa) {
  grid <- 0:min(min(ua), M)
  if (joint_kappa) {
    m <- kappa + mu
    grid <- grid[(m - grid) > 0 & (m - grid) <= M]
    lp <- vapply(grid, function(k)
      sum(na * dnbinom(ua - k, size = r, mu = m - k, log = TRUE)), numeric(1))
  } else {
    lp <- vapply(grid, function(k)
      sum(na * dnbinom(ua - k, size = r, mu = mu, log = TRUE)), numeric(1))
  }
  k_new <- if (length(grid) == 1) grid else
    grid[sample.int(length(grid), 1, prob = exp(lp - max(lp)))]
  list(kappa = k_new, mu = if (joint_kappa) m - k_new else mu)
}

# Random walk step on rtilde, with lambda integrated out ---------------------
draw_rtilde <- function(ya, na, mu, rtilde, sd_rw) {
  cand <- rnorm(1, rtilde, sd_rw)
  if (cand <= 0 || cand >= 0.5) return(list(rtilde = rtilde, acc = 0))
  ll <- function(rt) sum(na * dnbinom(ya, size = (1 - rt) / rt, mu = mu, log = TRUE))
  if (log(runif(1)) < ll(cand) - ll(rtilde)) list(rtilde = cand, acc = 1)
  else list(rtilde = rtilde, acc = 0)
}

# lambda, then mu given lambda: 1/mu is Gamma(r n - 1, rate r * sum(lambda)) truncated to [1/M, Inf)
draw_mu <- function(ya, na, mu, r, M) {
  n <- sum(na)
  sum_lambda <- rgamma(1, shape = sum(na * ya) + n * r, scale = mu / (mu + r))
  a <- n * r - 1
  b <- r * sum_lambda
  upper <- pgamma(1 / M, a, rate = b, lower.tail = FALSE)
  if (upper < 1e-12) return(M)
  1 / qgamma(runif(1) * upper, a, rate = b, lower.tail = FALSE)
}

# Random walk on log(alpha), with the weights integrated out -----------------
draw_alpha <- function(alpha, nj, step = 1, shape0 = 1, rate0 = 200) {
  J <- length(nj); n <- sum(nj)
  lt <- function(a) (shape0 - 1) * log(a) - rate0 * a +
    lgamma(J * a) - lgamma(n + J * a) + sum(lgamma(nj + a) - lgamma(a))
  cand <- exp(log(alpha) + step * rnorm(1))
  # log(cand) - log(alpha) is the Jacobian of the transformation to log(alpha)
  if (log(runif(1)) < lt(cand) - lt(alpha) + log(cand) - log(alpha)) list(alpha = cand, acc = 1)
  else list(alpha = alpha, acc = 0)
}

# The sampler ----------------------------------------------------------------
snb_sampler <- function(y, J, burnin = 5000, nkeep = 10000, sparse = (J > 1),
                        alpha = if (sparse) 0.005 else 1, joint_kappa = TRUE) {
  tab <- table(y); u <- as.numeric(names(tab)); cnt <- as.numeric(tab)
  m_u <- length(u); M <- max(y)

  kappa <- rep(0, J)
  rtilde <- rep(1 / 11, J)                    # r = 10
  mu <- rep(max(mean(y), 1), J)
  w <- rep(1 / J, J)
  sd_rw <- rep(0.1, J)
  nuj <- matrix(cnt, m_u, 1)                  # counts by unique value and component
  if (J > 1) nuj <- t(sapply(cnt, function(cc) rmultinom(1, cc, rep(1 / J, J))))

  store <- list(mu = matrix(NA, nkeep, J), rtilde = matrix(NA, nkeep, J),
                kappa = matrix(NA, nkeep, J), w = matrix(NA, nkeep, J),
                nj = matrix(NA, nkeep, J), alpha = rep(NA, nkeep))
  batch_acc <- batch_n <- keep_acc <- keep_n <- rep(0, J)
  n_batch <- 0; alpha_acc <- 0

  for (it in 1:(burnin + nkeep)) {
    # Block 1: allocation, with lambda integrated out
    if (J > 1) {
      lp <- matrix(-Inf, m_u, J)
      for (j in 1:J) {
        id <- u >= kappa[j]
        lp[id, j] <- log(w[j]) + dnbinom(u[id] - kappa[j], size = (1 - rtilde[j]) / rtilde[j],
                                         mu = mu[j], log = TRUE)
      }
      pr <- exp(lp - apply(lp, 1, max))
      for (i in 1:m_u) nuj[i, ] <- rmultinom(1, cnt[i], pr[i, ])
    }
    nj <- colSums(nuj)

    for (j in 1:J) {
      if (nj[j] == 0) {                       # an empty component is drawn from its priors
        kappa[j] <- sample.int(M + 1, 1) - 1
        mu[j] <- M * (1 - runif(1))
        rtilde[j] <- 0.5 * (1 - runif(1))
        next
      }
      act <- nuj[, j] > 0; ua <- u[act]; na <- nuj[act, j]
      # Block 2: kappa
      s <- draw_kappa(ua, na, kappa[j], mu[j], (1 - rtilde[j]) / rtilde[j], M, joint_kappa)
      kappa[j] <- s$kappa; mu[j] <- s$mu
      ya <- ua - kappa[j]
      # Block 3: rtilde
      s <- draw_rtilde(ya, na, mu[j], rtilde[j], sd_rw[j])
      rtilde[j] <- s$rtilde
      batch_acc[j] <- batch_acc[j] + s$acc; batch_n[j] <- batch_n[j] + 1
      if (it > burnin) { keep_acc[j] <- keep_acc[j] + s$acc; keep_n[j] <- keep_n[j] + 1 }
      # Blocks 4 and 5: lambda and mu
      mu[j] <- draw_mu(ya, na, mu[j], (1 - rtilde[j]) / rtilde[j], M)
    }

    # Block 6: alpha
    if (sparse) {
      s <- draw_alpha(alpha, nj); alpha <- s$alpha
      if (it > burnin) alpha_acc <- alpha_acc + s$acc
    }
    # Block 7: weights
    g <- rgamma(J, nj + alpha, 1)
    if (!all(is.finite(g)) || sum(g) == 0) g <- nj + 1e-10
    w <- g / sum(g)

    # Batch adaptation of the candidate sd for rtilde, in burn-in only (target 0.3 to 0.4)
    if (it <= burnin && it %% 50 == 0) {
      n_batch <- n_batch + 1; dl <- min(0.1, 1 / sqrt(n_batch))
      for (j in 1:J) if (batch_n[j] > 0) {
        rate <- batch_acc[j] / batch_n[j]
        if (rate > 0.4) sd_rw[j] <- sd_rw[j] * exp(dl)
        else if (rate < 0.3) sd_rw[j] <- sd_rw[j] * exp(-dl)
      }
      batch_acc[] <- 0; batch_n[] <- 0
    }

    if (it > burnin) {
      i <- it - burnin
      store$mu[i, ] <- mu; store$rtilde[i, ] <- rtilde; store$kappa[i, ] <- kappa
      store$w[i, ] <- w; store$nj[i, ] <- nj; store$alpha[i] <- alpha
    }
  }
  store$acc <- keep_acc / pmax(keep_n, 1)
  store$alpha_acc <- alpha_acc / nkeep
  store
}

# Number and location of the modes for one draw, as in Section 3.2 -----------
# Non-empty components only, weights renormalized, every integer up to just above the largest
# component mode. A plateau counts as one mode, and all of its locations are returned.
snb_modes <- function(mu, rtilde, kappa, w, nj) {
  act <- nj > 0
  mu <- mu[act]; rtilde <- rtilde[act]; kappa <- kappa[act]; w <- w[act] / sum(w[act])
  r <- (1 - rtilde) / rtilde
  top <- max(kappa + floor((r - 1) * mu / r)) + 1
  yg <- 0:top
  pm <- rep(0, length(yg))
  for (j in seq_along(r)) {
    id <- yg >= kappa[j]
    pm[id] <- pm[id] + w[j] * dnbinom(yg[id] - kappa[j], size = r[j], mu = mu[j])
  }
  s <- sign(diff(c(0, pm)))
  nz <- which(s != 0)
  L <- 0; locs <- integer(0)
  for (t in seq_along(nz)) {
    i <- nz[t]
    if (s[i] == 1) {
      if (t == length(nz)) nxt <- length(s) + 1
      else if (s[nz[t + 1]] == -1) nxt <- nz[t + 1]
      else next
      L <- L + 1
      locs <- c(locs, yg[i:(nxt - 1)])
    }
  }
  list(L = L, locs = locs)
}

# Summaries ------------------------------------------------------------------
summarize_modes <- function(o) {
  K <- nrow(o$mu)
  res <- lapply(1:K, function(i) snb_modes(o$mu[i, ], o$rtilde[i, ], o$kappa[i, ], o$w[i, ], o$nj[i, ]))
  L <- vapply(res, function(x) x$L, numeric(1))
  locs <- unlist(lapply(res, function(x) x$locs))
  list(number = round(table(L) / K, 3), location = round(table(locs) / K, 3))
}

# Parameters for draws with the most frequent number of non-empty components, sorted by kappa then mu.
# Labels are only needed for this table: the modes do not depend on them.
summarize_parameters <- function(o) {
  K0 <- rowSums(o$nj > 0); Khat <- as.numeric(names(which.max(table(K0))))
  keep <- which(K0 == Khat)
  arr <- lapply(keep, function(i) {
    act <- which(o$nj[i, ] > 0)
    act <- act[order(o$kappa[i, act], o$mu[i, act])]
    c(o$mu[i, act], (1 - o$rtilde[i, act]) / o$rtilde[i, act], o$kappa[i, act], o$w[i, act] / sum(o$w[i, act]))
  })
  arr <- do.call(rbind, arr)
  nm <- c(paste0("mu_", 1:Khat), paste0("r_", 1:Khat), paste0("kappa_", 1:Khat), paste0("weight_", 1:Khat))
  tab <- t(apply(arr, 2, function(v) c(mean = mean(v), quantile(v, c(0.025, 0.5, 0.975)))))
  rownames(tab) <- nm
  list(Khat = Khat, K0 = round(table(K0) / nrow(o$nj), 3), table = round(tab, 3))
}

simulate_snb <- function(n, w, kappa, mu, r) {
  z <- sample.int(length(w), n, replace = TRUE, prob = w)
  kappa[z] + rnbinom(n, size = r[z], mu = mu[z])
}

# Recovery check on designs D1 to D4 -----------------------------------------
run_check <- function(joint_kappa = TRUE, burnin = 5000, nkeep = 10000, chains = 2) {
  designs <- list(
    D1 = list(n = 2000,  w = 1,          kappa = 5,       mu = 10,      r = 3,       J = 1, sparse = FALSE),
    D2 = list(n = 10000, w = c(.5, .5),  kappa = c(0, 50), mu = c(15, 1), r = c(10, 2), J = 2, sparse = FALSE),
    D3 = list(n = 2000,  w = c(.6, .4),  kappa = c(0, 10), mu = c(2, 5),  r = c(5, 3),  J = 2, sparse = FALSE),
    D4 = list(n = 2000,  w = c(.6, .4),  kappa = c(0, 10), mu = c(2, 5),  r = c(5, 3),  J = 6, sparse = TRUE))
  for (d in seq_along(designs)) {
    D <- designs[[d]]
    set.seed(3000 + d)
    y <- simulate_snb(D$n, D$w, D$kappa, D$mu, D$r)
    truth <- snb_modes(D$mu, 1 / (1 + D$r), D$kappa, D$w, rep(1, length(D$w)))
    cat(sprintf("\n===== %s, n = %d, J = %d, joint_kappa = %s. True modes: %s\n", names(designs)[d], D$n, D$J,
                joint_kappa, paste(truth$locs, collapse = ", ")))
    cat("True values: mu =", D$mu, "| r =", D$r, "| kappa =", D$kappa, "| weight =", D$w, "\n")
    for (ch in 1:chains) {
      set.seed(7000 + 10 * d + ch)
      t0 <- Sys.time()
      o <- snb_sampler(y, D$J, burnin, nkeep, sparse = D$sparse, alpha = if (D$sparse) 0.005 else 1, joint_kappa = joint_kappa)
      tm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("Chain %d (%.0f s), rtilde acceptance: %s\n", ch, tm, paste(round(o$acc[o$acc > 0], 2), collapse = " ")))
      p <- summarize_parameters(o); m <- summarize_modes(o)
      cat("  Non-empty components:", paste(names(p$K0), p$K0, sep = ": ", collapse = "  "), "\n")
      cat("  Number of modes:", paste(names(m$number), m$number, sep = ": ", collapse = "  "), "\n")
      keepl <- m$location[m$location > 0.01]
      cat("  Pr(value is a mode):", paste(names(keepl), keepl, sep = ": ", collapse = "  "), "\n")
      print(p$table)
    }
  }
}

# Runs only with Rscript, not when the file is sourced
if (sys.nframe() == 0) run_check(joint_kappa = TRUE)
