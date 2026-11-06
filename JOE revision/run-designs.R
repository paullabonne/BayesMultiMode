designs <- list(
  D1 = list(n = 2000L, K = 1L, weight = 1, shift = 5, mu = 10, shape = 3,
            data_seed = 202601L, chain_seed = 31001L),
  D2 = list(n = 10000L, K = 2L, weight = c(0.5, 0.5), shift = c(0, 50),
            mu = c(15, 1), shape = c(10, 2), data_seed = 202602L, chain_seed = 31002L),
  D3 = list(n = 2000L, K = 2L, weight = c(0.6, 0.4), shift = c(0, 10),
            mu = c(2, 5), shape = c(5, 3), data_seed = 202603L, chain_seed = 31003L),
  D4 = list(n = 2000L, K = 6L, weight = c(0.6, 0.4), shift = c(0, 10),
            mu = c(2, 5), shape = c(5, 3), data_seed = 202603L, chain_seed = 31004L)
)

true_pmf <- function(design, grid) {
  Reduce(`+`, lapply(seq_along(design$weight), function(component) {
    design$weight[component] * dnbinom(grid - design$shift[component],
                                     size = design$shape[component], mu = design$mu[component])
  }))
}

draw_pmf <- function(draw, count, grid) {
  Reduce(`+`, lapply(seq_len(count), function(component) {
    shape <- 1 / draw[[paste0("rtilde", component)]] - 1
    draw[[paste0("eta", component)]] * dnbinom(
      grid - draw[[paste0("kappa", component)]], size = shape,
      mu = draw[[paste0("mu", component)]])
  }))
}

pmf_modes <- function(grid, probability) {
  grid[probability > 0 & probability > c(-Inf, head(probability, -1)) &
         probability >= c(tail(probability, -1), -Inf)]
}

aligned_parameters <- function(draw, count, true_count) {
  weight <- draw[paste0("eta", seq_len(count))]
  shift <- draw[paste0("kappa", seq_len(count))]
  mu <- draw[paste0("mu", seq_len(count))]
  rtilde <- draw[paste0("rtilde", seq_len(count))]
  shape <- 1 / rtilde - 1
  selected <- order(weight, decreasing = TRUE)[seq_len(true_count)]
  selected <- selected[order((shift + mu)[selected], shift[selected])]
  extras <- setdiff(order(weight, decreasing = TRUE), selected)
  cbind(weight, shift, mu, shape, rtilde, p = shape / (shape + mu),
        component_mean = shift + mu)[c(selected, extras), , drop = FALSE]
}

parameter_table <- function(fit, design) {
  parameter_names <- c("weight", "shift", "mu", "shape", "rtilde", "p", "component_mean")
  true_count <- length(design$weight)
  aligned <- vapply(seq_len(nrow(fit$mcmc)), function(draw) {
    as.vector(aligned_parameters(fit$mcmc[draw, ], fit$K, true_count))
  }, numeric(fit$K * length(parameter_names)))
  truth <- cbind(design$weight, design$shift, design$mu, design$shape,
                 1 / (1 + design$shape), design$shape / (design$shape + design$mu),
                 design$shift + design$mu)
  if (fit$K > true_count) {
    truth <- rbind(truth, cbind(rep(0, fit$K - true_count),
                               matrix(NA_real_, fit$K - true_count, 6)))
  }
  intervals <- t(apply(aligned, 1, quantile, probs = c(0.025, 0.975)))
  result <- data.frame(component = rep(seq_len(fit$K), length(parameter_names)),
                       parameter = rep(parameter_names, each = fit$K),
                       truth = as.vector(truth), estimate = rowMeans(aligned),
                       lower = intervals[, 1], upper = intervals[, 2])
  result[order(result$component, match(result$parameter, parameter_names)), ]
}

run_design <- function(name, nb_iter = 10000L, burnin = 5000L) {
  design <- designs[[name]]
  set.seed(design$data_seed)
  component <- sample.int(length(design$weight), design$n, replace = TRUE, prob = design$weight)
  observations <- design$shift[component] + rnbinom(design$n,
    size = design$shape[component], mu = design$mu[component])
  set.seed(design$chain_seed)
  started <- proc.time()[["elapsed"]]
  fit <- if (design$K == 1L) {
    priors <- BayesMultiMode:::check_priors(list(), "shifted_neg_binomial", observations)
    draws <- BayesMultiMode:::gibbs_SFM_shifted_neg_binomial(observations, K = design$K,
      nb_iter = nb_iter, nb_burnin = burnin, priors = priors, print = TRUE)
    parameters <- draws[, colnames(draws) != "loglik", drop = FALSE]
    list(data = observations, K = design$K, mcmc_all = parameters,
         mcmc = parameters[seq.int(burnin + 1L, nb_iter), , drop = FALSE],
         loglik = draws[, "loglik"])
  } else {
    BayesMultiMode::bayes_fit(observations, K = design$K, dist = "shifted_neg_binomial",
                              nb_iter = nb_iter, burnin = burnin, print = TRUE)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  stopifnot(nrow(fit$mcmc) == nb_iter - burnin, nrow(fit$mcmc_all) == nb_iter,
            all(is.finite(fit$mcmc)), all(is.finite(fit$loglik)))
  weights <- fit$mcmc[, paste0("eta", seq_len(fit$K)), drop = FALSE]
  stopifnot(all(weights >= 0), max(abs(rowSums(weights) - 1)) < 1e-10)
  grid <- seq.int(0, max(200, observations, design$shift +
    qnbinom(1 - 1e-10, size = design$shape, mu = design$mu)))
  probability <- vapply(seq_len(nrow(fit$mcmc)), function(draw) {
    draw_pmf(fit$mcmc[draw, ], fit$K, grid)
  }, numeric(length(grid)))
  stopifnot(all(is.finite(probability)), all(probability >= 0),
            all(colSums(probability) <= 1 + 1e-10))
  truth <- true_pmf(design, grid)
  estimated <- rowMeans(probability)
  stopifnot(abs(sum(truth) - 1) < 1e-8)
  selected <- unique(round(seq(1, nrow(fit$mcmc), length.out = 10)))
  likelihood <- vapply(selected, function(draw) {
    sum(log(draw_pmf(fit$mcmc[draw, ], fit$K, observations)))
  }, numeric(1))
  likelihood_error <- max(abs(likelihood - tail(fit$loglik, nrow(fit$mcmc))[selected]))
  stopifnot(likelihood_error < 1e-6)
  parameter_draws_mean <- vapply(seq_len(nrow(fit$mcmc)), function(draw) {
    parameters <- aligned_parameters(fit$mcmc[draw, ], fit$K, length(design$weight))
    sum(parameters[, "weight"] * parameters[, "component_mean"])
  }, numeric(1))
  list(name = name, design = design, fit = fit, nb_iter = nb_iter, burnin = burnin,
       elapsed = elapsed, grid = grid, truth = truth, estimated = estimated,
       lower = apply(probability, 1, quantile, probs = 0.025),
       upper = apply(probability, 1, quantile, probs = 0.975),
       parameters = parameter_table(fit, design),
       tv = (sum(abs(estimated - truth)) + abs(sum(estimated - truth))) / 2,
       cdf_error = max(abs(cumsum(estimated - truth))),
       tail_mass = max(0, 1 - sum(estimated)), likelihood_error = likelihood_error,
       mean_ess = posterior::ess_bulk(matrix(parameter_draws_mean, ncol = 1)),
       loglik_ess = posterior::ess_bulk(matrix(tail(fit$loglik, nrow(fit$mcmc)), ncol = 1)),
       extra_weight = if (fit$K > length(design$weight)) {
         apply(weights, 1, function(weight) sum(sort(weight, decreasing = TRUE)[
           seq.int(length(design$weight) + 1, fit$K)]))
       } else rep(0, nrow(weights)),
       session = sessionInfo())
}

plot_design <- function(result, path) {
  grDevices::png(path, width = 1500, height = 850, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.6, 4.8, 4.3, 1), mgp = c(2.8, 0.8, 0), las = 1)
  grid <- result$grid
  observed <- tabulate(result$fit$data + 1L, nbins = length(grid)) / result$design$n
  xmax <- max(result$fit$data, grid[which(cumsum(result$truth) >= 0.9995)[1]],
              grid[which(cumsum(result$estimated) >= min(0.9995, sum(result$estimated)))[1]])
  visible <- grid <= xmax
  graphics::plot(grid, result$truth, type = "n", xlim = c(0, xmax),
                 ylim = c(0, 1.15 * max(observed[visible], result$upper[visible], result$truth[visible])),
                 xlab = "Count", ylab = "Probability mass",
                 main = sprintf("%s: true and estimated distribution", result$name))
  graphics::mtext(sprintf("n = %s | K = %d | one chain | %s retained draws",
    format(result$design$n, big.mark = ","), result$design$K,
    format(nrow(result$fit$mcmc), big.mark = ",")), side = 3, line = 0.4, cex = 0.85)
  graphics::abline(h = graphics::axTicks(2), col = "#EEEEEE")
  graphics::segments(grid, 0, grid, observed, col = "#CDD1D4", lwd = 4)
  graphics::polygon(c(grid, rev(grid)), c(result$lower, rev(result$upper)),
                    col = grDevices::adjustcolor("#0072B2", alpha.f = 0.18), border = NA)
  graphics::lines(grid, result$truth, col = "#202020", lwd = 2)
  graphics::points(grid, result$truth, col = "#202020", pch = 16, cex = 0.42)
  graphics::lines(grid, result$estimated, col = "#0072B2", lwd = 2, lty = 2)
  graphics::legend("topright", c("True PMF", "Posterior mean PMF", "95% pointwise interval", "Simulated data"),
                   col = c("#202020", "#0072B2", "#BDD9E9", "#CDD1D4"),
                   lty = c(1, 2, NA, NA), lwd = c(2, 2, NA, NA),
                   pch = c(16, NA, 15, 15), pt.cex = c(0.6, 1, 1.5, 1.5), bty = "n", cex = 0.85)
}

format_estimate <- function(value) {
  if (is.na(value)) return("N/A")
  if (value != 0 && abs(value) < 0.0001) return(sprintf("%.3e", value))
  sprintf("%.4f", value)
}

write_report <- function(results, directory) {
  lines <- c("# D1-D4: single-chain simulation results", "",
    "One simulated dataset and one MCMC chain per design. D3 and D4 use the same observations.", "",
    "## Method", "",
    "Each chain uses 10,000 iterations: 5,000 burn-in and 5,000 retained draws, without thinning. Fits use the current local shifted negative binomial sampler and its default priors. No chain was selected from multiple runs.", "",
    "D1 calls the underlying sampler directly because the public mixture-object constructor rejects `K = 1`; D2-D4 use `bayes_fit`. All four use the same sampling algorithm and default prior settings (`a0 = 1`, `A0 = 200`, initial MH scale `sig = 0.1`).", "",
    "The estimated distribution is the mean of the mixture PMF over all retained draws, including all six fitted components in D4. Shading gives 95% pointwise posterior intervals for the PMF, not predictive intervals for individual observations. Lines connect integer counts for readability.", "",
    "Parameter estimates are posterior means, with equal-tailed 95% credible intervals. Shape and success probability are transformed in each draw before averaging: `r = 1/rtilde - 1` and `p = r/(r + mu)`. Here `mu` is the unshifted negative binomial mean; the component mean is `kappa + mu`.", "",
    "For component comparisons, retain the largest-weight components (one in D1, two otherwise) and order them by `kappa + mu` in each draw. D4 extras are ranked by decreasing weight in each draw. This is a reporting convention, not evidence that exactly two components were recovered. Extras have true weight zero and undefined true shape/location parameters.", "",
    "One chain and one dataset per design do not establish convergence or repeated-sampling recovery. No between-chain R-hat is available. Single-chain bulk ESS values below describe only the sampled trajectory.", "",
    "## Distribution comparison", "",
    "TV pools the probability above the evaluation grid into one tail bin; it is a lower bound on full-support total variation. CDF error is the maximum on the grid. Modes below belong to the true and posterior-mean PMFs, not the posterior distribution of mode counts.", "",
    "| Case | n | K | True modes | Estimated modes | TV (tail pooled) | Max CDF error | Mean ESS | Log-likelihood ESS | Runtime (s) |",
    "|---|---:|---:|---|---|---:|---:|---:|---:|---:|")
  for (result in results) {
    lines <- c(lines, sprintf("| %s | %d | %d | %s | %s | %.4f | %.4f | %.0f | %.0f | %.1f |",
      result$name, result$design$n, result$design$K,
      paste(pmf_modes(result$grid, result$truth), collapse = ", "),
      paste(pmf_modes(result$grid, result$estimated), collapse = ", "),
      result$tv, result$cdf_error, result$mean_ess, result$loglik_ess, result$elapsed))
  }
  worst <- results[[which.max(vapply(results, function(result) result$tv, numeric(1)))]]
  lines <- c(lines, "", sprintf("%s has the largest distribution error among these runs (TV %.4f). The parameter tables below also show why a close distribution fit does not establish parameter recovery.",
    worst$name, worst$tv))
  labels <- c(weight = "Weight pi", shift = "Shift kappa", mu = "NB mean mu", shape = "Shape r",
              rtilde = "rtilde", p = "Success probability p", component_mean = "Component mean kappa + mu")
  for (result in results) {
    image <- paste0("results-", result$name, ".png")
    plot_design(result, file.path(directory, image))
    lines <- c(lines, "", paste("##", result$name), "",
      sprintf("![%s true and estimated probability mass functions](%s)", result$name, image), "",
      sprintf("Data seed: `%d`. Chain seed: `%d`. PMF evaluation grid: 0-%d. Estimated mass beyond the grid: %.3g.",
        result$design$data_seed, result$design$chain_seed, max(result$grid), result$tail_mass))
    shifts <- subset(result$parameters, parameter == "shift" & component <= length(result$design$weight))
    missed <- shifts$component[shifts$truth < shifts$lower | shifts$truth > shifts$upper]
    if (length(missed)) {
      lines <- c(lines, "", sprintf("**Parameter recovery warning:** the true shift lies outside the sampled 95%% interval for component(s) %s. Narrow or degenerate intervals here do not establish accurate recovery or adequate exploration of the shift parameter.",
        paste(missed, collapse = ", ")))
    }
    if (is.finite(result$loglik_ess) && result$loglik_ess < 100) {
      lines <- c(lines, "", sprintf("**Mixing warning:** log-likelihood bulk ESS is only %.1f out of %d retained draws. Treat the reported means and intervals as descriptions of this run, not reliable converged posterior summaries.",
        result$loglik_ess, nrow(result$fit$mcmc)))
    }
    lines <- c(lines, "", "| Component | Parameter | True | Estimated mean | 95% credible interval |",
      "|---|---|---:|---:|---|")
    for (row in seq_len(nrow(result$parameters))) {
      entry <- result$parameters[row, ]
      label <- if (entry$component <= length(result$design$weight)) as.character(entry$component) else {
        paste("Extra", entry$component - length(result$design$weight))
      }
      lines <- c(lines, sprintf("| %s | %s | %s | %s | [%s, %s] |", label,
        labels[[entry$parameter]], format_estimate(entry$truth),
        format_estimate(entry$estimate), format_estimate(entry$lower), format_estimate(entry$upper)))
    }
    if (result$design$K > length(result$design$weight)) {
      interval <- quantile(result$extra_weight, c(0.025, 0.975))
      lines <- c(lines, "", sprintf("Combined weight of the four extras: posterior mean %s; 95%% credible interval [%s, %s]. Their individual parameter summaries describe weight-ranked slots, not persistent component identities.",
        format_estimate(mean(result$extra_weight)), format_estimate(interval[1]), format_estimate(interval[2])))
    }
  }
  lines <- c(lines, "", "## Reproduction and checks", "",
    "Run from the repository root:", "", "```r", "source(\"JOE revision/run-designs.R\")",
    "run_designs()", "```", "",
    "To regenerate only the report and plots from saved chains, without running any new chains:", "",
    "```r", "write_report(readRDS(\"JOE revision/results.rds\"), \"JOE revision\")", "```", "",
    "The runner saves observations, full chains, retained draws, summaries, and R session details in [results.rds](results.rds). Seeds and settings are recorded above. It checks finite draws, normalized weights, PMF non-negativity and mass, and agreement between stored and independently recomputed log-likelihoods at ten retained draws per case. These numerical checks do not certify convergence.", "",
    sprintf("Report generated: %s. %s. BayesMultiMode %s loaded from the local working tree, including uncommitted sampler changes.",
      format(Sys.Date()), R.version.string, as.character(utils::packageVersion("BayesMultiMode"))))
  writeLines(lines, file.path(directory, "results.md"))
}

run_designs <- function(cases = names(designs)) {
  stopifnot(length(cases) > 0, all(cases %in% names(designs)), !anyDuplicated(cases))
  pkgload::load_all(".", quiet = TRUE)
  directory <- "JOE revision"
  cache <- file.path(directory, "results.rds")
  results <- if (file.exists(cache)) readRDS(cache) else list()
  for (name in cases) {
    cat("\nRunning", name, "with one chain\n")
    results[[name]] <- run_design(name)
    results <- results[intersect(names(designs), names(results))]
    saveRDS(results, cache)
    write_report(results, directory)
    cat(sprintf("%s complete: TV = %.4f, CDF error = %.4f, elapsed = %.1f s\n",
      name, results[[name]]$tv, results[[name]]$cdf_error, results[[name]]$elapsed))
  }
  invisible(results)
}

if (sys.nframe() == 0L) {
  cases <- commandArgs(trailingOnly = TRUE)
  run_designs(if (length(cases)) cases else names(designs))
}