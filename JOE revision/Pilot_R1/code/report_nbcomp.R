# Report for run_nbcomp.R. Reference: our sampler with 50,000 + 50,000 iterations (res_full).
# IF = stored draws per effective draw (stored = every 5th of 10,000 kept iterations).
# ESS/s = effective draws per CPU second of the sampling loop (10,000 kept + 5,000 burn-in).
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
ALG <- c(ours = "Ours", mur_rw = "(mu, rtilde) RW", pr_rw = "(p, r) RW", zhou = "Zhou, his prior", zhou_mh = "Zhou-MH, our prior")
top_L <- function(L) as.numeric(names(which.max(table(L))))
rows <- list()
for (f in list.files(file.path(S, "res_nbcomp"), pattern = "rds$")) {
  X <- readRDS(file.path(S, "res_nbcomp", f))
  for (a in names(X)) {
    x <- X[[a]]; tk <- x$task; n <- if (tk$design == "default") 0 else 2000
    ref <- readRDS(file.path(S, "res_full", sprintf("%s_n%d_r%d_SNB_c%d.rds", tk$design, n, tk$rep, tk$chain)))
    rows[[length(rows) + 1]] <- data.frame(design = tk$design, rep = tk$rep, chain = tk$chain, alg = a, cpu = x$cpu,
      if_ll = x$if_[["loglik"]], if_L = x$if_[["modes"]], if_d = x$if_[["dispersion"]],
      eps_ll = x$ess_per_cpu_sec[["loglik"]], eps_d = x$ess_per_cpu_sec[["dispersion"]],
      gap = (mean(ref$loglik) - mean(x$loglik)) / sd(ref$loglik), same_top = top_L(x$L) == top_L(ref$L),
      pr_top = mean(x$L == top_L(ref$L)), pr_top_ref = mean(ref$L == top_L(ref$L)),
      acc1 = x$acc[1], acc2 = x$acc[2])
  }
}
A <- do.call(rbind, rows); saveRDS(A, file.path(S, "nbcomp_compare.rds"))
md <- function(v) median(v, na.rm = TRUE)
cat(format(Sys.time()), "| runs:", nrow(A), "\n")
cat("Medians over chains. IF: log-lik, modes, dispersion | ESS per CPU second: log-lik, dispersion | CPU s | far (gap > 1 sd) | same top number | Pr(ref top) [reference]\n")
for (d in c("G1", "G4", "G5", "G6", "default", "all")) {
  cat(d, "\n")
  for (a in names(ALG)) {
    z <- if (d == "all") A[A$alg == a, ] else A[A$design == d & A$alg == a, ]
    cat(sprintf("  %-19s | IF %5.1f %5.1f %6.1f | ESS/s %6.1f %6.2f | CPU %5.1f | far %d/%d | top %d/%d | Pr %.2f [%.2f] | acc %s\n",
                ALG[a], md(z$if_ll), md(z$if_L), md(z$if_d), md(z$eps_ll), md(z$eps_d), md(z$cpu), sum(z$gap > 1), nrow(z),
                sum(z$same_top), nrow(z), mean(z$pr_top), mean(z$pr_top_ref), paste(round(c(md(z$acc1), md(z$acc2)), 2), collapse = " ")))
  }
}
cat("\nRatio to ours of the median IF (log-lik, modes, dispersion) and of median ESS per CPU second (log-lik, dispersion), all chains\n")
for (a in names(ALG)[-1]) {
  r <- function(v, inv = FALSE) { x <- md(A[A$alg == a, v]) / md(A[A$alg == "ours", v]); x }
  cat(sprintf("  %-19s | IF x %4.1f %4.1f %5.1f | ESS/s x %4.2f %4.2f\n", ALG[a], r("if_ll"), r("if_L"), r("if_d"), r("eps_ll"), r("eps_d")))
}
