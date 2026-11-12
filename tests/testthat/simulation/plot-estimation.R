plot_simulation_method <- function(results, file) {
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  rows <- ceiling(length(results) / 2)
  panels <- matrix(seq_len(2 * rows), nrow = rows, byrow = TRUE)
  legend_panel <- 2 * rows + 1
  grDevices::png(file, width = 1800, height = 800 * rows + 100, res = 150)
  on.exit(grDevices::dev.off())
  graphics::layout(rbind(panels, rep(legend_panel, 2)), heights = c(rep(7, rows), 1))

  for (result in results) {
    experiment <- result$experiment
    grid <- experiment$grid
    truth <- simulation_density(experiment, grid)
    step <- diff(grid)[1]
    limits <- range(c(result$observations, grid[c(
      which(cumsum(truth) * step >= 0.001)[1],
      which(cumsum(truth) * step >= 0.999)[1]
    )]))
    breaks <- if (step == 1) {
      seq(-0.5, max(result$observations) + 0.5, by = 1)
    } else seq(min(result$observations), max(result$observations), length.out = 61)
    histogram <- graphics::hist(result$observations, breaks = breaks, plot = FALSE)
    visible <- grid >= limits[1] & grid <= limits[2]
    curves <- c(list(truth), Filter(Negate(is.null), result$densities))
    maximum <- max(c(histogram$density, unlist(lapply(curves, function(curve) curve[visible]))))
    graphics::par(mar = c(6.5, 4.2, 3.5, 1), mgp = c(2.4, 0.7, 0))
    graphics::plot(histogram, freq = FALSE, col = "#DFE3E5", border = "white",
                    xlim = limits, ylim = c(0, 1.08 * maximum), xlab = "Observation",
                    ylab = if (step == 1) "Probability" else "Density",
                    main = if (!is.null(experiment$title)) experiment$title else {
                      if (experiment$modes == 1L) "Unimodal" else "Bimodal"
                    })
    for (chain in seq_along(result$densities)) {
      density <- result$densities[[chain]]
      if (!is.null(density)) {
        graphics::lines(grid, density, col = colors[chain], lty = chain + 1, lwd = 1.8)
      }
    }
    graphics::lines(grid, truth, col = "#191919", lwd = 2.3)
    labels <- vapply(seq_along(result$densities), function(chain) {
      density <- result$densities[[chain]]
      if (is.null(density)) return(paste0("C", chain, ": no fit"))
      distance <- (sum(abs(density - truth)) * step + abs(sum(density - truth) * step)) / 2
      sprintf("C%d TV=%.3f", chain, distance)
    }, character(1))
    graphics::mtext(paste(labels, collapse = "   "), side = 1, line = 4, cex = 0.7)
    failed <- which(nzchar(result$errors))
    if (length(failed)) {
      errors <- unique(result$errors[failed])
      message <- paste(strwrap(paste(errors, collapse = "; "), width = 85), collapse = "\n")
      graphics::mtext(message, side = 1, line = 5.2, cex = 0.58, col = "#A32F20")
    }
    fit_settings <- utils::modifyList(list(K = 10L, count = 600L), result$fit)
    graphics::mtext(paste(experiment$dist, "| K =", fit_settings$K, "| n =", fit_settings$count),
                      side = 3, line = 0.4,
                      cex = 0.8, col = "#555555")
  }

  if (length(results) %% 2) {
    graphics::par(mar = rep(0, 4))
    graphics::plot.new()
  }
  graphics::par(mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend("center", legend = c("Simulated data", "Truth", paste("Chain", 1:4)),
                    col = c("#B5BEC3", "#191919", colors), pch = c(15, rep(NA, 5)),
                    lty = c(0, 1, 2:5), lwd = c(0, 2.3, rep(1.8, 4)),
                    bty = "n", horiz = TRUE, cex = 0.95)
}