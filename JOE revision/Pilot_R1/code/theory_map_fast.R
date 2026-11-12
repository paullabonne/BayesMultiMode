# Exact smallest TV between one shifted NB(kappa, mu, r) and any unshifted NB mixture (mixed Poisson),
# with the closed-form log-convexity condition kappa (r - 1) <= 2.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
src <- readLines(file.path(S, "theory_lp_fast.R")); src <- src[1:(grep("^logconvex_viol", src) - 1)]
eval(parse(text = src))
out <- expand.grid(mu = c(1, 5, 20), r = c(1.5, 2, 3, 5, 10, 50), kappa = c(1, 2, 3, 5, 10, 20))
out$cond <- out$kappa * (out$r - 1) <= 2
out$lower <- out$upper <- out$crude <- NA
for (i in seq_len(nrow(out))) {
  p <- snb_pmf(1, out$kappa[i], out$mu[i], out$r[i]); m <- min_tv(p, step = 0.5)
  out$lower[i] <- m$lower; out$upper[i] <- m$upper; out$crude[i] <- tv_floor(p); cat(i, out$kappa[i], out$mu[i], out$r[i], round(m$lower, 3), "
"); saveRDS(out, file.path(S, "theory_map.rds"))
}
saveRDS(out, file.path(S, "theory_map.rds"))
for (mu in unique(out$mu)) {
  cat(sprintf("\nmu = %g: exact min TV (rows kappa, columns r); * marks kappa (r - 1) <= 2\n", mu))
  x <- out[out$mu == mu, ]
  cat(sprintf("%8s", "kappa"), sprintf("%8s", paste0("r=", unique(x$r))), "\n")
  for (k in unique(x$kappa)) {
    z <- x[x$kappa == k, ]
    cat(sprintf("%8d", k), sprintf("%7.3f%s", z$lower, ifelse(z$cond, "*", " ")), "\n")
  }
}
cat("\nmax gap between attained and certified lower:", round(max(out$upper - out$lower), 4), "\n")
