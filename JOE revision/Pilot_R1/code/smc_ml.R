# Table 3 pilot: log-predictive likelihood = sum over observations of log p(y_t | y_1, ..., y_{t-1})
# (Geweke, 2005, Section 2.6.2), i.e. the log marginal likelihood, estimated by sequential Monte Carlo
# that adds one observation at a time (Chopin, 2002). The incremental weights are the one-step-ahead
# predictive likelihoods. Particles start from the prior of Section 2.2; after each resampling they are
# moved by `sweeps` sweeps of the sampler of NB_Comparisons.R ("ours" or "sp") on the data seen so far.
# Usage: Rscript smc_ml.R <data> <model> <N> <sweeps> <seed>
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R"))
a <- commandArgs(TRUE); dname <- a[1]; model <- a[2]; N <- as.integer(a[3]); sweeps <- as.integer(a[4]); seed <- as.integer(a[5])
sp <- model == "sp"; J <- 10
y <- real_data(dname); n <- length(y); M <- max(y); u <- sort(unique(y)); m_u <- length(u)
dir.create(file.path(S, "res_smc"), showWarnings = FALSE)
set.seed(seed); y <- y[sample.int(n)]                       # random order; the sum is order-invariant in theory
iu <- match(y, u)

# Particles: matrices N x J for mu, rtilde, kappa, w; vector alpha
draw_prior <- function(N) {
  al <- rgamma(N, 1, rate = 200)
  lw <- t(sapply(al, function(s) { g <- log(rgamma(J, s + 1)) + log(runif(J)) / s; g - (max(g) + log(sum(exp(g - max(g))))) }))
  list(alpha = al, w = exp(lw), kappa = matrix(sample.int(M + 1, N * J, TRUE) - 1, N, J),
       mu = matrix(M * (1 - runif(N * J)), N, J), rtilde = matrix(if (sp) 0 else 0.5 * (1 - runif(N * J)), N, J))
}
pred <- function(P, yv) {                                     # p(yv | theta_k) for every particle
  out <- numeric(N)
  for (j in 1:J) { ok <- P$kappa[, j] <= yv
    out[ok] <- out[ok] + P$w[ok, j] * dnbinom(yv - P$kappa[ok, j], size = (1 - P$rtilde[ok, j]) / P$rtilde[ok, j], mu = P$mu[ok, j]) }
  out
}
move <- function(th, cnt) {                                   # `sweeps` sweeps of the sampler for one particle
  mu <- th$mu; rtilde <- th$rtilde; kappa <- th$kappa; w <- th$w; alpha <- th$alpha
  idu <- which(cnt > 0); uu <- u[idu]; cc <- cnt[idu]; mm <- length(uu)
  for (s in 1:sweeps) {
    lp <- matrix(-Inf, mm, J)
    for (j in 1:J) { id <- uu >= kappa[j]; if (any(id) && w[j] > 0)
      lp[id, j] <- log(w[j]) + dnbinom(uu[id] - kappa[j], size = (1 - rtilde[j]) / rtilde[j], mu = mu[j], log = TRUE) }
    pr <- exp(lp - apply(lp, 1, max)); nuj <- matrix(0, mm, J)
    for (i in 1:mm) nuj[i, ] <- rmultinom(1, cc[i], pr[i, ])
    nj <- colSums(nuj)
    for (j in 1:J) {
      if (nj[j] == 0) { kappa[j] <- sample.int(M + 1, 1) - 1; mu[j] <- M * (1 - runif(1)); rtilde[j] <- if (sp) 0 else 0.5 * (1 - runif(1)); next }
      act <- nuj[, j] > 0; ua <- uu[act]; na <- nuj[act, j]; r <- (1 - rtilde[j]) / rtilde[j]
      k <- draw_kappa(ua, na, kappa[j], mu[j], r, M); kappa[j] <- k$kappa; mu[j] <- k$mu
      ya <- ua - kappa[j]
      if (sp) mu[j] <- draw_mu_pois(ya, na, M) else {
        rtilde[j] <- draw_rtilde(ya, na, mu[j], rtilde[j], 0.1)$value
        mu[j] <- draw_mu(ya, na, mu[j], (1 - rtilde[j]) / rtilde[j], M) }
    }
    alpha <- draw_alpha(alpha, nj)
    g <- rgamma(J, nj + alpha, 1); if (!all(is.finite(g)) || sum(g) == 0) g <- nj + 1e-10
    w <- g / sum(g)
  }
  list(mu = mu, rtilde = rtilde, kappa = kappa, w = w, alpha = alpha)
}

t0 <- proc.time()
P <- draw_prior(N); logW <- rep(0, N); lml <- 0; cnt <- rep(0, m_u); nmove <- 0; lp_t <- numeric(n)
for (t in 1:n) {
  l <- pred(P, y[t])
  W <- exp(logW - max(logW)); W <- W / sum(W)
  lp_t[t] <- log(sum(W * l)); lml <- lml + lp_t[t]
  logW <- logW + log(pmax(l, 1e-300)); cnt[iu[t]] <- cnt[iu[t]] + 1
  W <- exp(logW - max(logW)); ess <- 1 / sum((W / sum(W))^2)
  if (ess < N / 2 || t == n) {
    pos <- (runif(1) + 0:(N - 1)) / N; id <- findInterval(pos, cumsum(W / sum(W))) + 1; id[id > N] <- N
    for (k in 1:N) {
      th <- move(list(mu = P$mu[id[k], ], rtilde = P$rtilde[id[k], ], kappa = P$kappa[id[k], ], w = P$w[id[k], ], alpha = P$alpha[id[k]]), cnt)
      P2 <- if (k == 1) list(mu = matrix(0, N, J), rtilde = matrix(0, N, J), kappa = matrix(0, N, J), w = matrix(0, N, J), alpha = numeric(N)) else P2
      P2$mu[k, ] <- th$mu; P2$rtilde[k, ] <- th$rtilde; P2$kappa[k, ] <- th$kappa; P2$w[k, ] <- th$w; P2$alpha[k] <- th$alpha
    }
    P <- P2; logW <- rep(0, N); nmove <- nmove + 1
  }
}
tt <- proc.time() - t0
cat(sprintf("%s %s N=%d sweeps=%d seed=%d | n=%d | log-predictive likelihood %.1f | moves %d | CPU %.0f s\n",
            dname, model, N, sweeps, seed, n, lml, nmove, tt["user.self"] + tt["sys.self"]))
saveRDS(list(data = dname, model = model, N = N, sweeps = sweeps, seed = seed, lml = lml, lp_t = lp_t, y = y, nmove = nmove,
             cpu = unname(tt["user.self"] + tt["sys.self"])),
        file.path(S, "res_smc", sprintf("%s_%s_N%d_s%d_seed%d.rds", dname, model, N, sweeps, seed)))
