# Smallest TV distance between a target pmf and any mixed Poisson distribution.
# Every unshifted NB mixture (any J) is mixed Poisson, so this is a lower bound for all of them,
# and it is attained in the limit because Gamma mixtures are dense among mixing distributions.
# Upper value: TV of a feasible mixing distribution on a lambda grid (quadprog, tiny ridge).
# Certified lower bound (LP duality): for any g in [0,1]^Y,
#   min TV >= sum_y g_y p_y - sup_lambda sum_y g_y Pois(y; lambda), sup over a fine grid.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R"))
library(quadprog)
min_tv <- function(p, step = 0.25) {
  keep <- which(p > 1e-6); ys <- 0:(max(keep) - 1 + 10); pp <- p[ys + 1]
  lam <- seq(0, 1.3 * max(ys), by = step)
  A <- outer(ys, lam, function(y, l) dpois(y, l))
  K <- length(lam); Y <- length(ys); eps <- 1e-7
  Amat <- cbind(c(rep(1, K), rep(0, Y)),                 # sum w = 1
                rbind(t(A), diag(Y)),                    # A w + u >= p
                diag(K + Y))                             # w, u >= 0
  bvec <- c(1, pp, rep(0, K + Y))
  s <- solve.QP(Dmat = diag(eps, K + Y), dvec = -c(rep(0, K), rep(1, Y)), Amat = Amat, bvec = bvec, meq = 1)
  w <- pmax(s$solution[1:K], 0); w <- w / sum(w)
  fit <- as.vector(A %*% w)
  upper <- sum(pmax(pp - fit, 0)) + sum(pmax(p[-(ys + 1)], 0))   # TV of this feasible mixture
  g <- pmin(pmax(s$Lagrangian[2:(Y + 1)], 0), 1)
  fine <- seq(0, 1.3 * max(ys) + 20, by = 0.01)
  lower <- sum(g * pp) - max(vapply(fine, function(l) sum(g * dpois(ys, l)), 0))
  list(upper = upper, lower = lower, fit = fit, ys = ys, w = w, lam = lam)
}
logconvex_viol <- function(p, ys = 1:60) {
  # mixed Poisson: q(y) = y! p(y) is log-convex, i.e. p(y)^2 <= p(y-1) p(y+1) (y+1)/y
  r <- sapply(ys, function(y) p[y + 1]^2 / (p[y] * p[y + 2] * (y + 1) / y))
  ys[which(r > 1.0001)]
}
cat("design | crude limit | min TV over unshifted mixtures: certified lower, attained | y where y!p(y) is not log-convex\n")
res <- list()
for (d in names(DESIGNS)) {
  p <- DESIGNS[[d]]$pmf; m <- min_tv(p); res[[d]] <- m
  cat(sprintf("%s | %.3f | %.3f  %.3f | %s\n", d, tv_floor(p), m$lower, m$upper, paste(head(logconvex_viol(p), 15), collapse = " ")))
}
saveRDS(res, file.path(S, "theory_lp.rds"))
cat("\nPoint mass at k: p(k) <= e^-k k^k / k!, so min TV = 1 - that\n")
for (k in c(1, 2, 5, 10, 20, 50)) cat(sprintf("  k = %d: bound %.3f, min TV %.3f\n", k, exp(-k) * k^k / factorial(k), 1 - exp(-k) * k^k / factorial(k)))
cat("\nOne NB with mean 2 and r = 10, shifted by k\n")
for (k in c(0, 2, 5, 10, 20)) {
  p <- snb_pmf(1, k, 2, 10); m <- min_tv(p)
  cat(sprintf("  k = %d: min TV lower %.3f attained %.3f | not log-convex at %s\n", k, m$lower, m$upper, paste(head(logconvex_viol(p), 8), collapse = " ")))
}
