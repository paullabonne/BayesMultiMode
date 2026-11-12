# Practical validation V1-V4 for the existing shifted negative binomial sampler.
# Protocol: JOE revision/experiment/02-bayesian-validation-experiments.md
# Run from the repository root: Rscript "JOE revision/experiment/run-validation.R"

pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(posterior))

out_dir <- file.path("JOE revision", "experiment", "validation-results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nb_iter <- 5000L
burnin <- 2500L
n_chains <- 4L
time_limit <- 55 * 60
run_start <- proc.time()[["elapsed"]]

experiments <- list(
  V1 = list(K = 1L, n = 1000L, weight = 1, shift = 5, mu = 8, shape = 3,
            data_seeds = c(20260611L, 20260612L)),
  V2 = list(K = 2L, n = 1000L, weight = c(0.6, 0.4), shift = c(0, 12), mu = c(2, 3),
            shape = c(5, 3), data_seeds = c(20260621L, 20260622L)),
  V3 = list(K = 2L, n = 1000L, weight = c(0.5, 0.5), shift = c(0, 1), mu = c(5, 5),
            shape = c(4, 4), data_seeds = c(20260631L, 20260632L)),
  V4 = list(K = 6L, n = 1000L, weight = c(0.6, 0.4), shift = c(0, 12), mu = c(2, 3),
            shape = c(5, 3), data_seeds = c(20260621L, 20260622L))  # same observations as V2
)

simulate_counts <- function(design, seed) {
  set.seed(seed)
  component <- sample.int(length(design$weight), design$n, replace = TRUE, prob = design$weight)
  design$shift[component] + rnbinom(design$n, size = design$shape[component], mu = design$mu[component])
}

true_pmf <- function(design, x) {
  Reduce(`+`, Map(function(w, s, m, r) w * dnbinom(x - s, size = r, mu = m),
                  design$weight, design$shift, design$mu, design$shape))
}

parameters_of <- function(mcmc, K) {
  list(eta = mcmc[, paste0("eta", seq_len(K)), drop = FALSE],
       kappa = mcmc[, paste0("kappa", seq_len(K)), drop = FALSE],
       mu = mcmc[, paste0("mu", seq_len(K)), drop = FALSE],
       shape = 1 / mcmc[, paste0("rtilde", seq_len(K)), drop = FALSE] - 1)
}

draw_pmf <- function(pars, draw, x) {
  Reduce(`+`, lapply(seq_len(ncol(pars$eta)), function(k)
    pars$eta[draw, k] * dnbinom(x - pars$kappa[draw, k], size = pars$shape[draw, k],
                                mu = pars$mu[draw, k])))
}

pmf_matrix <- function(pars, x) {
  vapply(seq_len(nrow(pars$eta)), function(draw) draw_pmf(pars, draw, x), numeric(length(x)))
}

# Modes on an integer grid; flat runs count once, boundaries are checked against zero padding.
pmf_modes <- function(x, p, tol = 1e-10) {
  padded <- c(0, p, 0)
  locations <- c(x[1] - 1, x, x[length(x)] + 1)
  step <- abs(diff(padded)) > tol * pmax(padded[-1], padded[-length(padded)])
  runs <- split(seq_along(padded), cumsum(c(TRUE, step)))
  value <- vapply(runs, function(i) padded[i[1]], numeric(1))
  n_run <- length(value)
  is_mode <- value > 0 & value > c(0, value[-n_run]) & value > c(value[-1], 0)
  unname(lapply(runs[is_mode], function(indices) as.integer(locations[indices])))
}

start_state <- function(y, K, seed) {
  set.seed(seed)
  km <- kmeans(y, centers = K, nstart = 30)
  centers <- as.vector(km$centers)
  centers[centers == 0] <- max(y) / 2
  list(mu = centers, size = as.vector(table(factor(km$cluster, levels = seq_len(K)))))
}

fit_chain <- function(y, K, seed) {
  start <- start_state(y, K, seed)
  set.seed(seed)
  began <- proc.time()[["elapsed"]]
  if (K == 1L) {
    priors <- check_priors(list(), "shifted_neg_binomial", y)
    draws <- gibbs_SFM_shifted_neg_binomial(y, K = K, nb_iter = nb_iter, nb_burnin = burnin,
                                            priors = priors, print = FALSE)
    loglik_full <- draws[, "loglik"]
    mcmc <- draws[(burnin + 1L):nb_iter, colnames(draws) != "loglik", drop = FALSE]
  } else {
    fit <- bayes_fit(y, K = K, dist = "shifted_neg_binomial", nb_iter = nb_iter,
                     burnin = burnin, print = FALSE)
    loglik_full <- fit$loglik
    mcmc <- fit$mcmc
  }
  elapsed <- proc.time()[["elapsed"]] - began
  stopifnot(nrow(mcmc) == nb_iter - burnin)
  loglik <- tail(loglik_full, nrow(mcmc))
  pars <- parameters_of(mcmc, K)
  rtilde <- mcmc[, paste0("rtilde", seq_len(K)), drop = FALSE]
  idx <- unique(round(seq(1, nrow(mcmc), length.out = 10)))
  recomputed <- vapply(idx, function(d) sum(log(draw_pmf(pars, d, y))), numeric(1))
  checks <- c(weight_sum_error = max(abs(rowSums(pars$eta) - 1)),
              min_weight = min(pars$eta),
              min_shift = min(pars$kappa),
              min_rtilde = min(rtilde),
              max_rtilde = max(rtilde),
              loglik_finite = as.numeric(all(is.finite(loglik_full))),
              loglik_recompute_error = max(abs(recomputed - loglik[idx])))
  checks <- c(checks, checks_pass = as.numeric(
    checks[["weight_sum_error"]] < 1e-10 && checks[["min_weight"]] >= 0 &&
      checks[["min_rtilde"]] > 0 && checks[["max_rtilde"]] < 1 &&
      checks[["loglik_finite"]] == 1 && checks[["loglik_recompute_error"]] < 1e-6))
  list(seed = seed, K = K, mcmc = mcmc, loglik = loglik, loglik_full = loglik_full,
       start = start, elapsed = elapsed, checks = checks)
}

recovery_rows <- function(pars, design, name, index) {
  K <- design$K
  n <- nrow(pars$eta)
  comp_mean <- pars$kappa + pars$mu
  ord <- matrix(vapply(seq_len(n), function(d) order(comp_mean[d, ]), integer(K)), ncol = K, byrow = TRUE)
  sort_draws <- function(M) {
    matrix(vapply(seq_len(n), function(d) M[d, ord[d, ]], numeric(K)), ncol = K, byrow = TRUE)
  }
  params <- list(weight = sort_draws(pars$eta), shift = sort_draws(pars$kappa),
                 mu_unshifted = sort_draws(pars$mu), shape = sort_draws(pars$shape),
                 component_mean = sort_draws(comp_mean))
  true_order <- order(design$shift + design$mu)
  truth <- list(weight = design$weight, shift = design$shift, mu_unshifted = design$mu,
                shape = design$shape, component_mean = design$shift + design$mu)
  truth <- lapply(truth, function(v) v[true_order])
  do.call(rbind, lapply(names(params), function(p) {
    M <- params[[p]]
    do.call(rbind, lapply(seq_len(K), function(k) {
      q <- quantile(M[, k], c(0.025, 0.975))
      data.frame(experiment = name, dataset = index, component = k, parameter = p,
                 truth = truth[[p]][k], posterior_mean = mean(M[, k]),
                 lower_95 = unname(q[1]), upper_95 = unname(q[2]),
                 covered = truth[[p]][k] >= q[1] && truth[[p]][k] <= q[2],
                 stringsAsFactors = FALSE)
    }))
  }))
}

plot_trace <- function(path, name, index, chains, K) {
  grDevices::png(path, width = 1500, height = 900, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.6, 3, 1), las = 1)
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  full <- lapply(chains, `[[`, "loglik_full")
  graphics::plot(NULL, xlim = c(1, length(full[[1]])), ylim = range(unlist(full)),
                 xlab = "Iteration", ylab = "Log likelihood",
                 main = sprintf("%s dataset %d: log likelihood (burn-in left of dashed line)", name, index))
  for (ch in seq_along(full)) graphics::lines(full[[ch]], col = cols[ch])
  graphics::abline(v = burnin, lty = 2)
  means <- lapply(chains, function(ch) {
    p <- parameters_of(ch$mcmc, K)
    rowSums(p$eta * (p$kappa + p$mu))
  })
  graphics::plot(NULL, xlim = c(1, length(means[[1]])), ylim = range(unlist(means)),
                 xlab = "Retained draw", ylab = "Mixture mean",
                 main = sprintf("%s dataset %d: mixture mean (retained draws)", name, index))
  for (ch in seq_along(means)) graphics::lines(means[[ch]], col = cols[ch])
  graphics::legend("topright", paste("chain", seq_along(chains)), col = cols[seq_along(chains)],
                   lty = 1, bty = "n", cex = 0.8)
}

plot_distribution <- function(path, name, index, y, grid, truth_p, est, band) {
  grDevices::png(path, width = 1500, height = 850, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.6, 4.8, 4.3, 1), las = 1)
  reach <- function(p) {
    i <- which(cumsum(p) >= 0.9995)
    if (length(i)) grid[i[1]] else max(grid)
  }
  xmax <- max(y, reach(truth_p), reach(est))
  visible <- grid <= xmax
  observed <- tabulate(y + 1L, nbins = length(grid)) / length(y)
  ymax <- 1.15 * max(observed[visible], band[2, visible], truth_p[visible])
  graphics::plot(grid, truth_p, type = "n", xlim = c(0, xmax), ylim = c(0, ymax),
                 xlab = "Count", ylab = "Probability mass",
                 main = sprintf("%s dataset %d: true and posterior-mean distribution", name, index))
  graphics::segments(grid, 0, grid, observed, col = "#CDD1D4", lwd = 4)
  graphics::polygon(c(grid, rev(grid)), c(band[1, ], rev(band[2, ])),
                    col = grDevices::adjustcolor("#0072B2", alpha.f = 0.18), border = NA)
  graphics::lines(grid, truth_p, col = "#202020", lwd = 2)
  graphics::lines(grid, est, col = "#0072B2", lwd = 2, lty = 2)
  graphics::legend("topright", c("True PMF", "Posterior mean PMF", "95% pointwise band", "Observed frequencies"),
                   col = c("#202020", "#0072B2", "#BDD9E9", "#CDD1D4"),
                   lty = c(1, 2, 1, 1), lwd = c(2, 2, 8, 4), bty = "n", cex = 0.85)
}

analyse_dataset <- function(name, design, index, y, chains) {
  K <- design$K
  mcmc <- do.call(rbind, lapply(chains, `[[`, "mcmc"))
  loglik <- unlist(lapply(chains, `[[`, "loglik"))
  per_chain <- nrow(chains[[1]]$mcmc)
  n_draws <- nrow(mcmc)
  chain_id <- rep(seq_along(chains), each = per_chain)
  pars <- parameters_of(mcmc, K)

  upper <- max(design$shift + qnbinom(1 - 1e-12, size = design$shape, mu = design$mu),
               pars$kappa + qnbinom(1 - 1e-12, size = pars$shape, mu = pars$mu), max(y))
  grid <- seq.int(0, ceiling(upper) + 5)
  P <- pmf_matrix(pars, grid)
  stopifnot(all(is.finite(P)), all(P >= 0))
  omitted <- 1 - colSums(P)
  excess <- max(colSums(P)) - 1

  truth_p <- true_pmf(design, grid)
  stopifnot(abs(sum(truth_p) - 1) < 1e-8)
  truth_peaks <- pmf_modes(grid, truth_p)
  truth_modes <- as.integer(unlist(truth_peaks))
  L0 <- length(truth_peaks)
  modes_list <- lapply(seq_len(n_draws), function(j) pmf_modes(grid, P[, j]))
  L_draw <- lengths(modes_list)
  L_ind <- as.numeric(L_draw == L0)
  est <- rowMeans(P)
  est_modes <- as.integer(unlist(pmf_modes(grid, est)))
  band <- apply(P, 1, quantile, probs = c(0.025, 0.975))
  tv <- 0.5 * (sum(abs(est - truth_p)) + abs((1 - sum(est)) - (1 - sum(truth_p))))
  cdf_error <- max(abs(cumsum(est - truth_p)))

  comp_mean <- pars$kappa + pars$mu
  mix_mean <- rowSums(pars$eta * comp_mean)
  mix_var <- rowSums(pars$eta * (pars$mu + pars$mu^2 / pars$shape + comp_mean^2)) - mix_mean^2
  x_sel <- sort(unique(c(0, truth_modes)))
  p_sel <- t(P[x_sel + 1, , drop = FALSE])
  quantities <- cbind(loglik, mix_mean, mix_var, p_sel, L_ind)
  colnames(quantities) <- c("loglik", "mixture_mean", "mixture_var",
                            paste0("pmf_at_", x_sel), "mode_count_is_true")

  Q <- ncol(quantities)
  arr <- array(NA_real_, dim = c(per_chain, length(chains), Q),
               dimnames = list(NULL, NULL, colnames(quantities)))
  for (q in seq_len(Q)) arr[, , q] <- matrix(quantities[, q], nrow = per_chain)
  summ <- summarise_draws(as_draws_array(arr), "rhat", "ess_bulk", "ess_tail")
  idx <- match(colnames(quantities), summ$variable)
  rhat <- summ$rhat[idx]
  ess_b <- summ$ess_bulk[idx]
  ess_t <- summ$ess_tail[idx]
  constant <- apply(arr, 3, function(v) all(v == v[1]))
  lag_acf <- function(v, lag) {
    if (all(v == v[1])) return(NA_real_)
    mean(vapply(seq_len(ncol(v)), function(ch) acf(v[, ch], lag.max = lag, plot = FALSE)$acf[lag + 1],
                numeric(1)))
  }
  ok <- (rhat < 1.01) & (ess_b > 400) & (ess_t > 400)
  ok[is.na(ok)] <- FALSE
  diagnostics <- data.frame(
    experiment = name, dataset = index, quantity = colnames(quantities),
    rhat = rhat, ess_bulk = ess_b, ess_tail = ess_t,
    lag1_acf = vapply(seq_len(Q), function(q) lag_acf(arr[, , q], 1), numeric(1)),
    lag10_acf = vapply(seq_len(Q), function(q) lag_acf(arr[, , q], 10), numeric(1)),
    status = ifelse(constant, "constant", ifelse(ok, "pass", "fail")),
    stringsAsFactors = FALSE)

  L_levels <- 0:max(L_draw)
  count_dist <- as.numeric(table(factor(L_draw, levels = L_levels))) / n_draws
  chain_mode_prob <- vapply(seq_along(chains), function(ch) mean(L_ind[chain_id == ch]), numeric(1))
  distribution <- data.frame(
    experiment = name, dataset = index, K = K,
    true_modes = paste(truth_modes, collapse = ";"), true_mode_count = L0,
    mean_pmf_modes = paste(est_modes, collapse = ";"),
    mode_count_distribution = paste(paste0(L_levels, "=", sprintf("%.3f", count_dist)), collapse = " "),
    prob_true_mode_count = mean(L_ind), tv = tv, max_cdf_error = cdf_error,
    mean_omitted_tail = mean(omitted), max_omitted_tail = max(omitted),
    max_pmf_sum_excess = excess, grid_max = max(grid), stringsAsFactors = FALSE)

  all_loc <- table(as.integer(unlist(modes_list)))
  cond_draws <- modes_list[L_draw == L0]
  cond_loc <- table(as.integer(unlist(cond_draws)))
  locs <- sort(unique(c(as.integer(names(all_loc)), as.integer(names(cond_loc)))))
  p_all <- as.numeric(all_loc[as.character(locs)]) / n_draws
  p_cond <- as.numeric(cond_loc[as.character(locs)]) / max(1, length(cond_draws))
  p_all[is.na(p_all)] <- 0
  p_cond[is.na(p_cond)] <- 0
  modes <- data.frame(experiment = name, dataset = index, location = locs,
                      prob_mode = p_all, prob_mode_given_true_count = p_cond,
                      stringsAsFactors = FALSE)

  chain_table <- do.call(rbind, lapply(seq_along(chains), function(ch) {
    ci <- chains[[ch]]
    data.frame(experiment = name, dataset = index, chain = ch, seed = ci$seed, K = K,
               nb_iter = nb_iter, burnin = burnin, runtime_sec = round(ci$elapsed, 1),
               start_mu = paste(round(ci$start$mu, 3), collapse = ";"),
               start_size = paste(ci$start$size, collapse = ";"),
               t(ci$checks), chain_mode_true_prob = chain_mode_prob[ch],
               stringsAsFactors = FALSE)
  }))

  plot_trace(file.path(out_dir, sprintf("trace-%s-%d.png", name, index)), name, index, chains, K)
  plot_distribution(file.path(out_dir, sprintf("distribution-%s-%d.png", name, index)),
                    name, index, y, grid, truth_p, est, band)

  recovery <- if (name %in% c("V1", "V2")) recovery_rows(pars, design, name, index) else NULL
  list(diagnostics = diagnostics, distribution = distribution, modes = modes,
       recovery = recovery, chains = chain_table)
}

results <- list()
for (e_index in seq_along(experiments)) {
  e_name <- names(experiments)[e_index]
  design <- experiments[[e_index]]
  for (d in seq_along(design$data_seeds)) {
    key <- paste0(e_name, ".", d)
    y <- simulate_counts(design, design$data_seeds[d])
    if (e_name == "V4") stopifnot(identical(y, simulate_counts(experiments$V2, experiments$V2$data_seeds[d])))
    chains <- list()
    for (ch in seq_len(n_chains)) {
      if (proc.time()[["elapsed"]] - run_start > time_limit) break
      seed <- 40000L + 100L * e_index + 10L * d + ch
      chains[[ch]] <- fit_chain(y, design$K, seed)
      message(sprintf("%s chain %d done (%.1f s)", key, ch, chains[[ch]]$elapsed))
    }
    if (length(chains) < n_chains) {
      message(key, ": time limit reached after ", length(chains), " chains; not analysed")
      next
    }
    results[[key]] <- analyse_dataset(e_name, design, d, y, chains)
    message(key, " analysed")
  }
}

diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
distribution <- do.call(rbind, lapply(results, `[[`, "distribution"))
modes <- do.call(rbind, lapply(results, `[[`, "modes"))
recovery <- do.call(rbind, lapply(results, `[[`, "recovery"))
chains_table <- do.call(rbind, lapply(results, `[[`, "chains"))

utils::write.csv(diagnostics, file.path(out_dir, "validation-diagnostics.csv"), row.names = FALSE)
utils::write.csv(distribution, file.path(out_dir, "validation-distribution.csv"), row.names = FALSE)
utils::write.csv(modes, file.path(out_dir, "validation-modes.csv"), row.names = FALSE)
utils::write.csv(recovery, file.path(out_dir, "validation-recovery.csv"), row.names = FALSE)
utils::write.csv(chains_table, file.path(out_dir, "validation-chains.csv"), row.names = FALSE)

cat("\nTotal elapsed (s):", round(proc.time()[["elapsed"]] - run_start, 1), "\n")
cat("\nChain checks passed:", sum(chains_table$checks_pass), "of", nrow(chains_table), "\n")
cat("\nDiagnostic status counts:\n"); print(table(diagnostics$status))
cat("\nQuantities not passing (fail):\n")
print(diagnostics[diagnostics$status == "fail", c("experiment", "dataset", "quantity", "rhat", "ess_bulk", "ess_tail")])
cat("\nDistribution summary:\n")
print(distribution[, c("experiment", "dataset", "true_modes", "mean_pmf_modes", "mode_count_distribution",
                       "prob_true_mode_count", "tv", "max_cdf_error", "max_omitted_tail")])
cat("\nRecovery:\n")
print(recovery)
