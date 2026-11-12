# Figure E.1: five samplers of NB_Comparisons.R (run_nbcomp.R, res_nbcomp), the run behind Table E.1.
# Writes figure_appE.pdf (for the manuscript) and figure_appE.png (to view).
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
A <- readRDS(file.path(S, "nbcomp_compare.rds"))
P <- c("ours", "mur_rw", "pr_rw", "zhou", "zhou_mh")
LAB <- c(ours = "Ours", mur_rw = expression((mu * "," ~ tilde(r)) ~ "random walks"), pr_rw = expression((p * "," ~ r) ~ "random walks"),
         zhou = "Zhou, own prior", zhou_mh = "Zhou as proposal, our prior")
COL <- c(ours = "#2a78d6", mur_rw = "#eb6834", pr_rw = "#1baf7a", zhou = "#eda100", zhou_mh = "#e87ba4")
PCH <- c(ours = 16, mur_rw = 17, pr_rw = 15, zhou = 18, zhou_mh = 25); LTY <- c(ours = 1, mur_rw = 2, pr_rw = 1, zhou = 3, zhou_mh = 4)
ink <- "#222222"; muted <- "#5F6670"; grid_col <- "#ECEEF0"
D <- c("G1", "G4", "G5", "G6", "default")
DL <- c("One negative\nbinomial", "Shifted second\ncomponent", "Narrow peak\nat 50", "Three close\npeaks", "Default\ndata")
md <- function(v) median(v, na.rm = TRUE)

draw <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(1, 1, 0.2))
  par(mar = c(3.8, 4.4, 2.8, 1), mgp = c(2.4, 0.6, 0), tcl = -0.25, las = 1, col.axis = muted, col.lab = muted, cex.axis = 0.9)
  ttl <- function(s, sub) { title(main = s, adj = 0, font.main = 1, cex.main = 1.05, col.main = ink, line = 1.1)
                            mtext(sub, side = 3, line = 0, adj = 0, cex = 0.8, col = muted) }
  # (a) trace of the dispersion summary, shifted second component, data set 1
  X <- readRDS(file.path(S, "res_nbcomp", "G4_r1_c1.rds")); sh <- c("ours", "pr_rw", "zhou")
  k <- 1:400; it <- 5000 + 5 * k
  tr <- lapply(sh, function(a) X[[a]]$disp[k]); names(tr) <- sh; rg <- range(unlist(tr))
  plot(NA, xlim = range(it), ylim = rg + c(-0.05, 0.05) * diff(rg), xlab = "Iteration", ylab = expression(Sigma[j] * pi[j] * log(r[j])), axes = FALSE)
  abline(h = pretty(rg), col = grid_col); axis(1, col = "#B8BDC3"); axis(2, col = NA, col.ticks = NA)
  for (a in sh) lines(it, tr[[a]], col = COL[a], lwd = 1.5, lty = LTY[a])
  ifs <- sapply(sh, function(a) X[[a]]$if_[["dispersion"]])
  ttl("(a) Dispersion summary, shifted second component",
      sprintf("Every 5th draw; inefficiency factor: ours %.1f, (p, r) %.0f, Zhou %.0f", ifs[1], ifs[2], ifs[3]))
  # (b)-(d) medians by design, log scale
  lolli <- function(v, lab, main, sub) {
    m <- sapply(P, function(a) sapply(D, function(d) md(A[A$design == d & A$alg == a, v])))
    yl <- range(m) * c(0.7, 1.6); base <- yl[1]
    plot(NA, xlim = c(0.5, 5.5), ylim = yl, log = "y", xlab = "", ylab = lab, axes = FALSE)
    at <- c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500); at <- at[at >= yl[1] & at <= yl[2]]
    abline(h = at, col = grid_col); axis(2, at = at, labels = format(at, drop0trailing = TRUE), col = NA, col.ticks = NA)
    axis(1, at = 1:5, labels = DL, col = "#B8BDC3", padj = 0.5, cex.axis = 0.82)
    off <- seq(-0.3, 0.3, length.out = 5)
    for (k in seq_along(P)) {
      segments(1:5 + off[k], base, 1:5 + off[k], m[, k], col = COL[P[k]], lwd = 2)
      points(1:5 + off[k], m[, k], pch = PCH[P[k]], col = COL[P[k]], bg = COL[P[k]], cex = 1.2)
    }
    ttl(main, sub)
  }
  lolli("if_ll", "Inefficiency factor (log scale)", "(b) Log-likelihood: inefficiency factor", "Median over chains; stored draws per effective draw")
  lolli("if_d", "Inefficiency factor (log scale)", "(c) Dispersion summary: inefficiency factor", "Median over chains; stored draws per effective draw")
  lolli("eps_d", "Effective draws per CPU second (log scale)", "(d) Dispersion summary: effective draws per CPU second", "Median over chains; higher is better")
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 5, bty = "n", legend = LAB[P], col = COL[P], pt.bg = COL[P], pch = PCH[P], lty = LTY[P], lwd = 2, text.col = ink, cex = 0.95, seg.len = 2.5)
}
cairo_pdf(file.path(S, "figure_appE.pdf"), width = 7, height = 6.2, pointsize = 10); draw(); dev.off()
png(file.path(S, "figure_appE.png"), width = 2600, height = 2120, res = 230, type = "cairo"); draw(); dev.off()
cat("done\n")
