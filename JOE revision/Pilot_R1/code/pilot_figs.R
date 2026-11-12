# Pilot versions of the application figures (Figures 5-7) for Referee 1, comment 5 and minor comments 2 and 6:
# larger fonts, thicker lines, no grid; number of modes on a common grid of at least three values for both models;
# mode locations marginally and given the two most probable numbers of modes. One PNG per data set.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
source(file.path(S, "s2_functions.R"))
rd <- function(d, m) lapply(1:4, function(c) readRDS(file.path(S, "res_pilot", sprintf("%s_%s_c%d.rds", d, m, c))))
ink <- "#222222"; muted <- "#5F6670"; bar <- "#C9CDD2"; c1 <- "#2a78d6"; c2 <- "#eb6834"; drawcol <- adjustcolor("#e34948", 0.08)
# Fixed symbol and colour for each number of modes, the same in every panel and figure:
# 1 square, 2 circle, 3 triangle, 4 diamond, 5 inverted triangle, 6 star, 7 cross (filled, coloured);
# 8 or more: open square, circle, triangle, diamond in black. Red is left for the posterior draws in the fit panel.
LPCH <- function(l) c(15, 16, 17, 18, 25, 8, 4, 0, 1, 2, 5)[pmin(l, 11)]
LCOL <- function(l) c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", rep(ink, 4))[pmin(l, 11)]
TT <- c(default ="Default data", inflation = "Inflation expectations (stand-in sample, March 2023)", dna = "DNA tandem repeats")
XL <- c(default = "Number of defaulted payments", inflation = "Expected inflation (%)", dna = "Number of repeat units")
for (d in names(TT)) {
  y <- real_data(d); R <- list(ours = rd(d, "ours"), sp = rd(d, "sp")); G <- R$ours[[1]]$G; xmax <- max(y)
  Ls <- lapply(R, function(ch) unlist(lapply(ch, `[[`, "L")))
  lg <- min(unlist(Ls)):max(unlist(Ls)); if (length(lg) < 3) lg <- seq(max(1, min(lg) - 1), length.out = 3)
  png(file.path(S, sprintf("pilot_fig_%s.png", d)), width = 3000, height = 1900, res = 230, type = "cairo")
  layout(matrix(c(1, 1, 1, 2:7), 3, 3, byrow = TRUE), heights = c(0.12, 1, 1))
  par(mar = c(0, 0, 0, 0)); plot.new(); text(0.01, 0.5, TT[d], adj = 0, cex = 1.6, col = ink)
  par(mar = c(4.2, 4.6, 2.6, 1), mgp = c(2.6, 0.7, 0), las = 1, tcl = -0.3, cex.axis = 1.15, cex.lab = 1.25, col.axis = ink, col.lab = ink)
  for (m in c("ours", "sp")) {
    ch <- R[[m]]; lab <- if (m == "ours") "Shifted negative binomial mixture" else "Shifted Poisson mixture"
    emp <- tabulate(y + 1, xmax + 1) / length(y)
    draws <- do.call(rbind, lapply(ch, `[[`, "pmf_draws")); pmean <- Reduce(`+`, lapply(ch, `[[`, "pmf")) / length(ch)
    yl <- c(0, max(emp, pmean[1:(xmax + 1)]) * 1.08)
    plot(NA, xlim = c(-0.5, xmax + 0.5), ylim = yl, xlab = XL[d], ylab = "Probability", bty = "l")
    rect(0:xmax - 0.4, 0, 0:xmax + 0.4, emp, col = bar, border = NA)
    for (k in seq(1, nrow(draws), by = 2)) lines(0:xmax, draws[k, 1:(xmax + 1)], col = drawcol, lwd = 1)
    lines(0:xmax, pmean[1:(xmax + 1)], col = ink, lwd = 2.2)
    title(main = paste0(lab, ": fit"), adj = 0, font.main = 1, cex.main = 1.25, col.main = ink)
    # number of modes, common grid
    pL <- sapply(lg, function(l) mean(Ls[[m]] == l))
    plot(NA, xlim = range(lg) + c(-0.6, 0.6), ylim = c(0, 1), xlab = "Number of modes", ylab = "Posterior probability", xaxt = "n", bty = "l")
    axis(1, at = lg); rect(lg - 0.35, 0, lg + 0.35, pL, col = LCOL(lg), border = NA)
    text(lg, pL, ifelse(pL >= .005, sprintf("%.2f", pL), ""), pos = 3, cex = 1.05, col = ink, xpd = NA)
    title(main = "Number of modes", adj = 0, font.main = 1, cex.main = 1.25, col.main = ink)
    # mode locations: marginal (bars) and given the two most probable numbers of modes (points)
    pm <- Reduce(`+`, lapply(ch, `[[`, "pr_mode")) / length(ch)
    tl <- as.numeric(names(sort(table(Ls[[m]]), decreasing = TRUE)))[1:min(2, length(unique(Ls[[m]])))]
    cond <- lapply(tl, function(l) Reduce(`+`, lapply(ch, function(s) { k <- match(as.character(l), colnames(s$pr_mode_L)); if (is.na(k)) 0 else s$pr_mode_L[, k] * sum(s$L == l) })) / sum(Ls[[m]] == l))
    plot(NA, xlim = c(-0.5, xmax + 0.5), ylim = c(0, 1.42), xlab = XL[d], ylab = "Pr(y is a mode)", bty = "l", yaxt = "n")
    axis(2, at = seq(0, 1, 0.2))
    rect(G - 0.4, 0, G + 0.4, pm, col = bar, border = NA)
    o <- order(tl); tl <- tl[o]; cond <- cond[o]                  # list in ascending number of modes
    for (k in seq_along(tl)) { w <- which(cond[[k]] >= 0.01); points(G[w], cond[[k]][w], col = LCOL(tl[k]), bg = LCOL(tl[k]), pch = LPCH(tl[k]), cex = 1.3) }
    legend("topleft", bty = "n", cex = 1.05, inset = c(0.01, 0), y.intersp = 0.95, legend = c("All draws", sprintf("Draws with %d mode%s (Pr %.2f)", tl, ifelse(tl == 1, "", "s"), sapply(tl, function(l) mean(Ls[[m]] == l)))),
           fill = c(bar, rep(NA, length(tl))), border = NA, pch = c(NA, LPCH(tl)), col = c(NA, LCOL(tl)), pt.bg = c(NA, LCOL(tl)), pt.cex = 1.3)
    title(main = "Mode locations", adj = 0, font.main = 1, cex.main = 1.25, col.main = ink)
  }
  dev.off()
}
cat("figures done\n")
