# Report for pilot_real.R (res_pilot) and smc_ml.R (res_smc).
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R"))
rd <- function(d, m, c) { f <- file.path(S, "res_pilot", sprintf("%s_%s_c%d.rds", d, m, c)); if (file.exists(f)) readRDS(f) else NULL }
top_L <- function(L) as.numeric(names(which.max(table(L))))
fmt <- function(p, G, cut = .1) { w <- which(p >= cut); paste(G[w], sprintf("%.2f", p[w]), sep = ":", collapse = "  ") }
for (d in c("default", "inflation", "dna")) for (m in c("ours", "sp")) {
  ch <- Filter(Negate(is.null), lapply(1:4, function(c) rd(d, m, c))); if (!length(ch)) next
  G <- ch[[1]]$G; L <- unlist(lapply(ch, `[[`, "L")); prL <- table(L) / length(L)
  pm <- Reduce(`+`, lapply(ch, `[[`, "pr_mode")) / length(ch)
  lppd <- mean(sapply(ch, `[[`, "lppd")); pw <- mean(sapply(ch, `[[`, "pwaic"))
  cat(sprintf("\n=== %s, %s (n = %d, %d chains) ===\n", d, if (m == "ours") "shifted NB" else "shifted Poisson", ch[[1]]$n, length(ch)))
  cat(sprintf(" WAIC %.0f (per chain %s) | split Rhat log-lik %.2f | chain tops %s\n", -2 * (lppd - pw),
              paste(round(sapply(ch, function(s) -2 * (s$lppd - s$pwaic))), collapse = " "), split_rhat(lapply(ch, `[[`, "loglik")),
              paste(sapply(ch, function(s) top_L(s$L)), collapse = " ")))
  cat(" Pr(number of modes):", paste(names(prL), sprintf("%.2f", prL), sep = ": ", collapse = "  "), "\n")
  cat(" Pr(y is a mode) >= .1:", fmt(pm, G), "\n")
  for (l in as.numeric(names(prL))[prL >= .05]) {
    pl <- Reduce(`+`, lapply(ch, function(s) { k <- match(as.character(l), colnames(s$pr_mode_L)); if (is.na(k)) 0 else s$pr_mode_L[, k] * sum(s$L == l) })) / sum(L == l)
    cat(sprintf("   given L = %g (Pr %.2f): %s\n", l, prL[as.character(l)], fmt(pl, G, .2)))
  }
  ifs <- function(v) median(sapply(ch, function(s) { x <- s[[v]]; length(x) / nb_ess(x) }), na.rm = TRUE)
  cat(sprintf(" Mixing (every 5th of 50,000; median over chains): IF log-lik %.1f | modes %.1f | dispersion %s | alpha %.1f\n",
              ifs("loglik"), ifs("L"), if (m == "sp") "n/a" else sprintf("%.1f", ifs("disp")), ifs("alpha")))
  cat(sprintf(" Acceptance: r tilde %s | alpha %.2f (step s = 1 on log alpha) | r tilde scales after burn-in, median %s, range %s\n",
              if (m == "sp") "n/a" else sprintf("%.2f", mean(sapply(ch, `[[`, "acc_rt"))), mean(sapply(ch, `[[`, "acc_alpha")),
              if (m == "sp") "n/a" else sprintf("%.3f", median(unlist(lapply(ch, `[[`, "sd_rt")))),
              if (m == "sp") "n/a" else paste(sprintf("%.3f", range(unlist(lapply(ch, `[[`, "sd_rt")))), collapse = "-")))
  cat(sprintf(" Non-empty components: median %g | alpha median %.4f | CPU s per chain %s\n", median(unlist(lapply(ch, `[[`, "K0"))),
              median(unlist(lapply(ch, `[[`, "alpha"))), paste(round(sapply(ch, `[[`, "cpu")), collapse = " ")))
}
cat("\n=== Table 3 pilot: log-predictive likelihood by SMC (N = 1,000 particles, 5 sweeps per move) ===\n")
fs <- list.files(file.path(S, "res_smc"), pattern = "N1000_s5", full.names = TRUE)
if (length(fs)) {
  X <- do.call(rbind, lapply(fs, function(f) { x <- readRDS(f); data.frame(data = x$data, model = x$model, seed = x$seed, lml = x$lml, moves = x$nmove, cpu = x$cpu) }))
  for (d in unique(X$data)) for (m in c("ours", "sp")) { z <- X[X$data == d & X$model == m, ]; if (!nrow(z)) next
    cat(sprintf(" %-9s %-4s | seeds %s | values %s | mean %.1f, range %.1f | moves %s | CPU s %s\n", d, m, paste(z$seed, collapse = ","),
                paste(sprintf("%.1f", z$lml), collapse = " "), mean(z$lml), diff(range(z$lml)), paste(z$moves, collapse = " "), paste(round(z$cpu), collapse = " "))) }
}
