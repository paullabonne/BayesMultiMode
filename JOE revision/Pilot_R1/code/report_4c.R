# Report for run_4c.R (R1 4c and 4d): mixture of shifted NB (ours) against mixture of shifted Poisson (sp).
# G1-G8: ours from res_full (same data sets, snb_s2_sampler.R), sp from res_4c. P1, P2: both from res_4c.
# Fit: TV to the true pmf, expected log score under the truth, WAIC. Modes: most probable number, Pr(true number).
# Efficiency: ESS of the log-likelihood and of the number of modes, from every 25th of 50,000 kept draws.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
src <- readLines(file.path(S, "run_4c.R"))
source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R"))
eval(parse(text = src[grep("^setup <- function", src):(grep("^run_4c_task <- function", src) - 1)])); setup()
top_L <- function(L) as.numeric(names(which.max(table(L))))
rd <- function(d, n, r, m, c) {
  f <- if (m == "SNB") { if (d %in% c("P1", "P2")) file.path(S, "res_4c", sprintf("%s_n%d_r%d_ours_c%d.rds", d, n, r, c)) else file.path(S, "res_full", sprintf("%s_n%d_r%d_SNB_c%d.rds", d, n, r, c)) }
       else file.path(S, "res_4c", sprintf("%s_n%d_r%d_sp_c%d.rds", d, n, r, c))
  readRDS(f)
}
one <- function(d, n, r, m, c = 1) {
  s <- rd(d, n, r, m, c); D <- DESIGNS[[d]]; Lt <- length(D$modes)
  secs <- if (!is.null(s$cpu)) s$cpu else s$secs
  data.frame(design = d, n = n, rep = r, model = m, chain = c, group = D$group,
             top = top_L(s$L), hit = top_L(s$L) == Lt, prL = mean(s$L == Lt), at = mean(s$at), near = mean(s$near),
             tv = 0.5 * sum(abs(s$pmf - D$pmf)), els = sum(D$pmf[D$pmf > 0] * log(pmax(s$pmf[D$pmf > 0], 1e-300))),
             waic = -2 * (s$lppd - s$pwaic), K0 = mean(s$K0),
             ess_ll = nb_ess(s$loglik), ess_L = nb_ess(s$L), secs = secs, timed = if (!is.null(s$cpu)) "cpu" else "elapsed")
}
des <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "P1", "P2")
A <- do.call(rbind, lapply(des, function(d) do.call(rbind, lapply(c(2000, 500), function(n) do.call(rbind, lapply(1:50, function(r) rbind(one(d, n, r, "SNB"), one(d, n, r, "SP"))))))))
saveRDS(A, file.path(S, "report_4c_chain1.rds"))

cat("TABLE A. Fit and modes, chain 1, 50 data sets per row (means; hit = share of data sets whose most probable number of modes is the true one)\n")
a <- aggregate(cbind(hit, prL, near, tv, els, K0) ~ design + n + model + group, A, mean)
a <- a[order(match(a$design, des), -a$n, a$model), ]
print(within(a, { hit <- round(hit, 2); prL <- round(prL, 2); near <- round(near, 2); tv <- round(tv, 3); els <- round(els, 4); K0 <- round(K0, 1) }), row.names = FALSE)

cat("\nTABLE B. Paired by data set, chain 1: SNB minus SP\n")
for (d in des) for (n in c(2000, 500)) {
  s <- A[A$design == d & A$n == n & A$model == "SNB", ]; p <- A[A$design == d & A$n == n & A$model == "SP", ]
  dw <- s$waic - p$waic; de <- s$els - p$els
  cat(sprintf(" %s n %4d | WAIC SNB lower in %2d/50, median diff %7.1f | expected log score SNB higher in %2d/50, mean diff %.4f | TV SNB %.3f SP %.3f | hit SNB %.2f SP %.2f\n",
              d, n, sum(dw < 0), median(dw), sum(de > 0), mean(de), mean(s$tv), mean(p$tv), mean(s$hit), mean(p$hit)))
}

cat("\nMost probable number of modes across the 50 data sets, n = 2000 and n = 500\n")
for (d in des) for (n in c(2000, 500)) for (m in c("SNB", "SP")) {
  x <- A[A$design == d & A$n == n & A$model == m, ]; t <- table(x$top)
  cat(sprintf(" %s n %4d %-3s true %d: %s\n", d, n, m, length(DESIGNS[[d]]$modes), paste(names(t), t, sep = " x", collapse = ", ")))
}

cat("\nWhere each model puts its modes, n = 2000: mean Pr(y is a mode), values at 0.2 or more\n")
for (d in des) for (m in c("SNB", "SP")) {
  p <- Reduce(`+`, lapply(1:50, function(r) rd(d, 2000, r, m, 1)$pr_mode)) / 50; w <- which(p >= 0.2)
  cat(sprintf(" %s %-3s true at %s: %s\n", d, m, paste(DESIGNS[[d]]$modes, collapse = ","), paste(GRID[w], round(p[w], 2), sep = ":", collapse = "  ")))
}

cat("\nTABLE C. Four chains, data sets 1 to 5\n")
for (d in des) for (n in c(2000, 500)) for (m in c("SNB", "SP")) {
  res <- t(sapply(1:5, function(r) {
    ch <- lapply(1:4, function(c) rd(d, n, r, m, c)); tops <- vapply(ch, function(s) top_L(s$L), 0)
    c(rhat = split_rhat(lapply(ch, function(s) s$loglik)), agree = length(unique(tops)) == 1,
      ll_range = diff(range(vapply(ch, function(s) mean(s$loglik), 0))))
  }))
  cat(sprintf(" %s n %4d %-3s | Rhat median %.2f max %.2f | chains agree on number of modes %d/5 | mean log-lik range max %.1f\n",
              d, n, m, median(res[, "rhat"]), max(res[, "rhat"]), sum(res[, "agree"]), max(res[, "ll_range"])))
}

cat("\nTABLE D. Efficiency, chain 1 (ESS from 2,000 summarised draws = every 25th of 50,000; IF = 50,000 / (25 x ESS) is a lower bound when ESS = 2,000)\n")
cat(" Timings: P1, P2 and all SP chains are CPU seconds from the same run; SNB on G1-G8 is elapsed seconds from the earlier res_full run.\n")
A$ess_ll_s <- A$ess_ll / A$secs; A$ess_L_s <- A$ess_L / A$secs
e <- aggregate(cbind(ess_ll, ess_L, secs, ess_ll_s, ess_L_s) ~ design + n + model, A, median, na.action = na.pass)
e <- e[order(match(e$design, des), -e$n, e$model), ]
print(within(e, { ess_ll <- round(ess_ll); ess_L <- round(ess_L); secs <- round(secs); ess_ll_s <- round(ess_ll_s, 1); ess_L_s <- round(ess_L_s, 1) }), row.names = FALSE)
