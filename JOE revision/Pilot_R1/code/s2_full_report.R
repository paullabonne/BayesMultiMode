# Report for the full R1 4a run. Reads whatever has finished; can be rerun at any time.
# Simulations: chain 1, only data sets where both models have finished, so the comparison is matched.
# Convergence: data sets 1 to 5, four chains.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R"))
fs <- list.files(file.path(S, "res_full"), pattern = "\\.rds$", full.names = TRUE)
R <- lapply(fs, readRDS)
tk <- do.call(rbind, lapply(R, function(s) data.frame(design = s$task$design, n = s$task$n, rep = s$task$rep, model = s$task$model, chain = s$task$chain)))
cat(length(R), "chains read,", format(Sys.time()), "| mean", round(mean(vapply(R, function(s) s$secs, 0))), "s per chain\n\n")

top_L <- function(L) as.numeric(names(which.max(table(L))))
one <- function(s, D) {
  Lt <- length(D$modes)
  data.frame(design = s$task$design, n = s$task$n, rep = s$task$rep, model = s$task$model,
             top = top_L(s$L), hit = top_L(s$L) == Lt, prL = mean(s$L == Lt), under = mean(s$L < Lt), over = mean(s$L > Lt),
             at = mean(s$at), near = mean(s$near), minnear = min(s$near),
             tv = 0.5 * sum(abs(s$pmf - D$pmf)), logscore = sum(D$pmf[D$pmf > 0] * log(pmax(s$pmf[D$pmf > 0], 1e-300))),
             K0 = mean(s$K0), secs = s$secs)
}
c1 <- which(tk$chain == 1 & tk$design %in% names(DESIGNS))
A <- do.call(rbind, lapply(c1, function(i) one(R[[i]], DESIGNS[[tk$design[i]]])))
key <- paste(A$design, A$n, A$rep)
A <- A[key %in% names(which(table(key) == 2)), ]
saveRDS(A, file.path(S, "full_chain1.rds"))

cat("TABLE 1. Chain 1, matched data sets\n")
a <- aggregate(cbind(hit, prL, under, over, at, near, minnear, tv, logscore, K0) ~ design + n + model, A, mean)
a$sets <- aggregate(rep ~ design + n + model, A, length)$rep
a$floor <- vapply(a$design, function(d) tv_floor(DESIGNS[[d]]$pmf), 0)
a <- a[order(a$design, -a$n, -xtfrm(a$model)), ]
print(within(a, { hit <- round(hit, 2); prL <- round(prL, 2); under <- round(under, 2); over <- round(over, 2); at <- round(at, 2); near <- round(near, 2); minnear <- round(minnear, 2); tv <- round(tv, 3); floor <- round(floor, 3); logscore <- round(logscore, 3); K0 <- round(K0, 1) })[, c("design", "n", "model", "sets", "hit", "prL", "under", "over", "at", "near", "minnear", "tv", "floor", "logscore", "K0")], row.names = FALSE)

cat("\nMost probable number of modes across data sets\n")
for (d in unique(A$design)) for (m in c("SNB", "NB")) {
  x <- A[A$design == d & A$model == m, ]; t <- table(x$top)
  cat(sprintf(" %s %-3s true %d: %s\n", d, m, length(DESIGNS[[d]]$modes), paste(names(t), t, sep = " x", collapse = ", ")))
}

cat("\nWhere each model puts its modes: mean Pr(y is a mode) across data sets, values at 0.2 or more\n")
for (d in unique(A$design)) for (m in c("SNB", "NB")) {
  id <- which(tk$design == d & tk$model == m & tk$chain == 1 & tk$n == 2000)
  p <- Reduce(`+`, lapply(R[id], function(s) s$pr_mode)) / length(id); w <- which(p >= 0.2)
  cat(sprintf(" %s %-3s true at %s: %s\n", d, m, paste(DESIGNS[[d]]$modes, collapse = ","), paste(GRID[w], round(p[w], 2), sep = ":", collapse = "  ")))
}

cat("\nTABLE 2. Pre-stated criteria, n = 2000\n")
for (d in unique(A$design)) {
  s <- a[a$design == d & a$n == 2000 & a$model == "SNB", ]; u <- a[a$design == d & a$n == 2000 & a$model == "NB", ]
  cat(sprintf(" %s (%s): right number SNB %.2f NB %.2f | TV SNB %.3f NB %.3f, ratio %.1f, gap %.3f | shift needed: %s | within .01: %s\n", d, DESIGNS[[d]]$group, s$hit, u$hit, s$tv, u$tv, u$tv / s$tv, u$tv - s$tv,
              (s$hit >= .9 && u$hit < .8) || (u$tv >= 2 * s$tv), abs(s$tv - u$tv) <= .01))
}

cat("\nTABLE 3. Four-chain checks, data sets 1 to 5\n")
for (d in names(DESIGNS)) for (n in sort(unique(A$n), decreasing = TRUE)) for (m in c("SNB", "NB")) {
  res <- lapply(1:5, function(r) {
    id <- which(tk$design == d & tk$n == n & tk$rep == r & tk$model == m)
    if (length(id) < 4) return(NULL)
    ch <- R[id]; tops <- vapply(ch, function(s) top_L(s$L), 0)
    prs <- vapply(ch, function(s) mean(s$L == length(DESIGNS[[d]]$modes)), 0)
    c(rhat = split_rhat(lapply(ch, function(s) s$loglik)), agree = length(unique(tops)) == 1, spread = max(prs) - min(prs),
      ll_range = diff(range(vapply(ch, function(s) mean(s$loglik), 0))))
  })
  res <- do.call(rbind, res); if (is.null(res)) next
  cat(sprintf(" %s n %d %-3s sets %d | Rhat median %.2f max %.2f, above 1.1: %d | all 4 chains same number of modes: %d | Pr(true number) spread max %.2f | mean log-lik range max %.1f\n", d, n, m, nrow(res), median(res[, "rhat"]), max(res[, "rhat"]), sum(res[, "rhat"] > 1.1), sum(res[, "agree"]), max(res[, "spread"]), max(res[, "ll_range"])))
}

cat("\nTABLE 4. Real data, four chains\n")
for (d in c("default", "dna", "inflation")) for (m in c("SNB", "NB")) {
  ch <- R[which(tk$design == d & tk$model == m)]; if (!length(ch)) next
  y <- real_data(d); L <- unlist(lapply(ch, function(s) s$L)); pm <- Reduce(`+`, lapply(ch, function(s) s$pr_mode)) / length(ch)
  lppd <- mean(vapply(ch, function(s) s$lppd, 0)); pw <- mean(vapply(ch, function(s) s$pwaic, 0))
  t <- round(table(L) / length(L), 2); t <- t[t >= .01]; w <- which(pm >= .1)
  cat(sprintf("%s %s n=%d chains %d: WAIC %.0f (per chain %s), Rhat %.2f, chain tops %s\n", d, m, length(y), length(ch), -2 * (lppd - pw), paste(round(vapply(ch, function(s) -2 * (s$lppd - s$pwaic), 0)), collapse = " "), split_rhat(lapply(ch, function(s) s$loglik)), paste(vapply(ch, function(s) top_L(s$L), 0), collapse = " ")))
  cat("   Pr(number of modes):", paste(names(t), t, sep = ": ", collapse = "  "), "\n")
  cat("   Pr(y is a mode):", paste(GRID[w], round(pm[w], 2), sep = ": ", collapse = "  "), "\n")
}
