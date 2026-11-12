# Comparison C1-C6: shifted NB vs unshifted NB vs shifted Poisson mixtures
# Protocol: JOE revision/experiment/03-shifted-nb-comparison-experiments.md
# Run from the repository root: Rscript "JOE revision/experiment/run-comparison.R"

pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(posterior))

default_out_dir <- file.path("JOE revision", "experiment", "comparison-results")
out_dir <- Sys.getenv("BAYESMULTIMODE_COMPARISON_OUT_DIR", default_out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_train <- as.numeric(Sys.getenv("BAYESMULTIMODE_COMPARISON_N_TRAIN", "500"))
if (length(n_train) != 1L || !is.finite(n_train) || n_train <= 0 || n_train != floor(n_train)) {
  stop("Training sample size must be a positive integer.")
}
if (n_train != 500L && identical(out_dir, default_out_dir)) {
  stop("A changed training sample size requires a separate output directory.")
}
n_test <- 5000L
K <- 6L
nb_iter <- 5000L
burnin <- 2500L
n_chains <- 4L
time_limit <- 55 * 60
run_start <- proc.time()[["elapsed"]]

models <- list(SNB = "shifted_neg_binomial", NB = "neg_binomial", SPois = "shifted_poisson")
selected_model <- Sys.getenv("BAYESMULTIMODE_COMPARISON_MODEL")
if (nzchar(selected_model) && !selected_model %in% names(models)) {
  stop("Unknown comparison model: ", selected_model)
}
spois_prior_name <- Sys.getenv("BAYESMULTIMODE_COMPARISON_SPOIS_PRIOR", "default")
spois_priors <- switch(spois_prior_name,
                      default = list(l0 = 5, L0 = 4),
                      broader = list(l0 = 1, L0 = 0.05),
                      stop("Unknown shifted-Poisson prior: ", spois_prior_name))

# NB(mu, r) uses mean mu and shape r; C6 NB(8,2) is read as mu = 8, r = 2.
settings <- list(
  C1 = list(kind = "snb", weight = 1, shift = 15, mu = 5, shape = 3),
  C2 = list(kind = "snb", weight = 1, shift = 0, mu = 20, shape = 3),
  C3 = list(kind = "snb", weight = c(0.5, 0.5), shift = c(10, 30), mu = c(5, 10), shape = c(100, 2)),
  C4 = list(kind = "pois", weight = c(0.5, 0.5), shift = c(10, 30), lambda = c(5, 10)),
  C5 = list(kind = "snb", weight = c(0.5, 0.5), shift = c(0, 1), mu = c(5, 5), shape = c(4, 4)),
  C6 = list(kind = "zi", p0 = 0.35, mu = 8, shape = 2)
)
selected_setting <- Sys.getenv("BAYESMULTIMODE_COMPARISON_SETTING")
if (nzchar(selected_setting) && !selected_setting %in% names(settings)) {
  stop("Unknown comparison setting: ", selected_setting)
}
setting_indices <- if (nzchar(selected_setting)) match(selected_setting, names(settings)) else seq_along(settings)
reuse_other_models <- nzchar(selected_model) && identical(out_dir, default_out_dir)

population_pmf <- function(s, x) {
  if (s$kind == "zi") {
    return(s$p0 * (x == 0) + (1 - s$p0) * dnbinom(x, size = s$shape, mu = s$mu))
  }
  out <- numeric(length(x))
  for (k in seq_along(s$weight)) {
    out <- out + s$weight[k] * if (s$kind == "pois") {
      dpois(x - s$shift[k], s$lambda[k])
    } else {
      dnbinom(x - s$shift[k], size = s$shape[k], mu = s$mu[k])
    }
  }
  out
}

pop_upper <- function(s) {
  if (s$kind == "zi") return(qnbinom(1 - 1e-12, size = s$shape, mu = s$mu))
  if (s$kind == "pois") return(max(s$shift + qpois(1 - 1e-12, s$lambda)))
  max(s$shift + qnbinom(1 - 1e-12, size = s$shape, mu = s$mu))
}

simulate_counts <- function(s, n, seed) {
  set.seed(seed)
  if (s$kind == "zi") {
    zero <- runif(n) < s$p0
    y <- rnbinom(n, size = s$shape, mu = s$mu)
    y[zero] <- 0L
    return(y)
  }
  component <- sample.int(length(s$weight), n, replace = TRUE, prob = s$weight)
  shift <- s$shift[component]
  if (s$kind == "pois") return(shift + rpois(n, s$lambda[component]))
  shift + rnbinom(n, size = s$shape[component], mu = s$mu[component])
}

parameters_of <- function(mcmc, dist, K) {
  idx <- seq_len(K)
  eta <- mcmc[, paste0("eta", idx), drop = FALSE]
  if (dist == "shifted_poisson") {
    return(list(type = "pois", eta = eta,
                kappa = mcmc[, paste0("kappa", idx), drop = FALSE],
                lambda = mcmc[, paste0("lambda", idx), drop = FALSE]))
  }
  kappa <- if (dist == "shifted_neg_binomial") {
    mcmc[, paste0("kappa", idx), drop = FALSE]
  } else {
    matrix(0, nrow(mcmc), K)
  }
  list(type = "nb", eta = eta, kappa = kappa,
       mu = mcmc[, paste0("mu", idx), drop = FALSE],
       shape = 1 / mcmc[, paste0("rtilde", idx), drop = FALSE] - 1)
}

draw_pmf <- function(pars, draw, x) {
  Reduce(`+`, lapply(seq_len(ncol(pars$eta)), function(k) {
    z <- x - pars$kappa[draw, k]
    pars$eta[draw, k] * if (pars$type == "pois") {
      dpois(z, pars$lambda[draw, k])
    } else {
      dnbinom(z, size = pars$shape[draw, k], mu = pars$mu[draw, k])
    }
  }))
}

pmf_matrix <- function(pars, x) {
  vapply(seq_len(nrow(pars$eta)), function(draw) draw_pmf(pars, draw, x), numeric(length(x)))
}

draw_upper <- function(pars) {
  if (pars$type == "pois") return(max(pars$kappa + qpois(1 - 1e-12, pars$lambda)))
  max(pars$kappa + qnbinom(1 - 1e-12, size = pars$shape, mu = pars$mu))
}

comp_moments <- function(pars) {
  if (pars$type == "pois") return(list(m = pars$kappa + pars$lambda, v = pars$lambda))
  list(m = pars$kappa + pars$mu, v = pars$mu + pars$mu^2 / pars$shape)
}

# Modes on an integer grid; flat runs count once, boundaries checked against zero padding.
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

fit_chain <- function(y, dist, seed) {
  start <- start_state(y, K, seed)
  set.seed(seed)
  began <- proc.time()[["elapsed"]]
  priors <- if (dist == "shifted_poisson") spois_priors else list()
  fit <- bayes_fit(y, K = K, dist = dist, nb_iter = nb_iter, burnin = burnin,
                   priors = priors, print = FALSE)
  elapsed <- proc.time()[["elapsed"]] - began
  mcmc <- fit$mcmc
  stopifnot(nrow(mcmc) == nb_iter - burnin)
  loglik_full <- fit$loglik
  loglik <- tail(loglik_full, nrow(mcmc))
  pars <- parameters_of(mcmc, dist, K)
  idx <- unique(round(seq(1, nrow(mcmc), length.out = 10)))
  recomputed <- vapply(idx, function(d) sum(log(draw_pmf(pars, d, y))), numeric(1))
  checks <- c(weight_sum_error = max(abs(rowSums(pars$eta) - 1)),
              min_weight = min(pars$eta),
              min_kappa = min(pars$kappa),
              loglik_finite = as.numeric(all(is.finite(loglik_full))),
              loglik_recompute_error = max(abs(recomputed - loglik[idx])))
  ok <- checks[["weight_sum_error"]] < 1e-10 && checks[["min_weight"]] >= 0 &&
    checks[["loglik_finite"]] == 1 && checks[["loglik_recompute_error"]] < 1e-6
  if (pars$type == "nb") {
    checks <- c(checks, min_shape = min(pars$shape), max_shape = max(pars$shape))
    ok <- ok && checks[["min_shape"]] > 0
  }
  checks <- c(checks, checks_pass = as.numeric(ok))
  list(seed = seed, mcmc = mcmc, loglik = loglik, loglik_full = loglik_full,
       start = start, elapsed = elapsed, checks = checks)
}

lag_acf <- function(v, lag) {
  if (all(v == v[1])) return(NA_real_)
  mean(vapply(seq_len(ncol(v)), function(ch) acf(v[, ch], lag.max = lag, plot = FALSE)$acf[lag + 1],
              numeric(1)))
}

chain_diagnostics <- function(quantities, per_chain, n_chains) {
  Q <- ncol(quantities)
  arr <- array(NA_real_, dim = c(per_chain, n_chains, Q),
               dimnames = list(NULL, NULL, colnames(quantities)))
  for (q in seq_len(Q)) arr[, , q] <- matrix(quantities[, q], nrow = per_chain)
  summ <- summarise_draws(as_draws_array(arr), "rhat", "ess_bulk", "ess_tail")
  idx <- match(colnames(quantities), summ$variable)
  rhat <- summ$rhat[idx]
  ess_b <- summ$ess_bulk[idx]
  ess_t <- summ$ess_tail[idx]
  constant <- apply(arr, 3, function(v) all(v == v[1]))
  ok <- (rhat < 1.01) & (ess_b > 400) & (ess_t > 400)
  ok[is.na(ok)] <- FALSE
  data.frame(quantity = colnames(quantities), rhat = rhat, ess_bulk = ess_b, ess_tail = ess_t,
             lag1_acf = vapply(seq_len(Q), function(q) lag_acf(arr[, , q], 1), numeric(1)),
             status = ifelse(constant, "constant", ifelse(ok, "pass", "fail")),
             stringsAsFactors = FALSE)
}

analyse_fit <- function(s_name, s, d, model_name, dist, y_train, y_test, chains) {
  mcmc <- do.call(rbind, lapply(chains, `[[`, "mcmc"))
  loglik <- unlist(lapply(chains, `[[`, "loglik"))
  per_chain <- nrow(chains[[1]]$mcmc)
  chain_id <- rep(seq_along(chains), each = per_chain)
  pars <- parameters_of(mcmc, dist, K)

  upper <- max(pop_upper(s), draw_upper(pars), max(y_test))
  grid <- seq.int(0, ceiling(upper) + 5)
  P <- pmf_matrix(pars, grid)
  stopifnot(all(is.finite(P)), all(P >= 0))
  omitted <- 1 - colSums(P)
  truth_p <- population_pmf(s, grid)
  stopifnot(abs(sum(truth_p) - 1) < 1e-8)

  pbar <- rowMeans(P)
  scores <- log(pbar[y_test + 1])
  chain_scores <- vapply(seq_along(chains), function(ch) {
    pc <- rowMeans(P[, chain_id == ch, drop = FALSE])
    mean(log(pc[y_test + 1]))
  }, numeric(1))
  tv <- 0.5 * (sum(abs(pbar - truth_p)) + abs((1 - sum(pbar)) - (1 - sum(truth_p))))
  q95 <- grid[which(cumsum(truth_p) >= 0.95)[1]]
  above <- grid > q95
  tail_error <- abs(sum(pbar[above]) - sum(truth_p[above]))

  truth_peaks <- pmf_modes(grid, truth_p)
  truth_modes <- as.integer(unlist(truth_peaks))
  L0 <- length(truth_peaks)
  modes_list <- lapply(seq_len(ncol(P)), function(j) pmf_modes(grid, P[, j]))
  L_draw <- lengths(modes_list)
  L_ind <- as.numeric(L_draw == L0)
  est_modes <- as.integer(unlist(pmf_modes(grid, pbar)))
  band <- apply(P, 1, quantile, probs = c(0.025, 0.975))

  mom <- comp_moments(pars)
  mix_mean <- rowSums(pars$eta * mom$m)
  mix_var <- rowSums(pars$eta * (mom$v + mom$m^2)) - mix_mean^2
  x_sel <- sort(unique(c(0, truth_modes)))
  quantities <- cbind(loglik, mix_mean, mix_var, t(P[x_sel + 1, , drop = FALSE]), L_ind)
  colnames(quantities) <- c("loglik", "mixture_mean", "mixture_var",
                            paste0("pmf_at_", x_sel), "mode_count_is_true")
  diag <- chain_diagnostics(quantities, per_chain, length(chains))

  L_levels <- 0:max(L_draw)
  count_dist <- as.numeric(table(factor(L_draw, levels = L_levels))) / ncol(P)
  chain_mode_prob <- vapply(seq_along(chains), function(ch) mean(L_ind[chain_id == ch]), numeric(1))

  all_loc <- table(as.integer(unlist(modes_list)))
  cond_draws <- modes_list[L_draw == L0]
  cond_loc <- table(as.integer(unlist(cond_draws)))
  locs <- sort(unique(c(as.integer(names(all_loc)), as.integer(names(cond_loc)))))
  p_all <- as.numeric(all_loc[as.character(locs)]) / ncol(P)
  p_cond <- as.numeric(cond_loc[as.character(locs)]) / max(1, length(cond_draws))
  p_all[is.na(p_all)] <- 0
  p_cond[is.na(p_cond)] <- 0

  summary_row <- data.frame(
    setting = s_name, dataset = d, model = model_name,
    n_train = length(y_train), n_test = length(y_test),
    tv = tv, mean_log_score = mean(scores), n_zero_density_obs = sum(!is.finite(scores)),
    chain_log_score_sd = sd(chain_scores), tail_error = tail_error,
    prob_true_mode_count = mean(L_ind), true_modes = paste(truth_modes, collapse = ";"),
    true_mode_count = L0, mean_pmf_modes = paste(est_modes, collapse = ";"),
    mode_count_distribution = paste(paste0(L_levels, "=", sprintf("%.3f", count_dist)), collapse = " "),
    max_omitted_tail = max(omitted), grid_max = max(grid),
    runtime_sec = round(sum(vapply(chains, `[[`, numeric(1), "elapsed")), 1),
    stringsAsFactors = FALSE)

  chain_table <- do.call(rbind, lapply(seq_along(chains), function(ch) {
    ci <- chains[[ch]]
    data.frame(setting = s_name, dataset = d, model = model_name, chain = ch, seed = ci$seed,
               n_train = length(y_train), n_test = length(y_test),
               K = K, nb_iter = nb_iter, burnin = burnin, runtime_sec = round(ci$elapsed, 1),
               prior_l0 = if (dist == "shifted_poisson") spois_priors$l0 else NA_real_,
               prior_L0 = if (dist == "shifted_poisson") spois_priors$L0 else NA_real_,
               start_mu = paste(round(ci$start$mu, 3), collapse = ";"),
               start_size = paste(ci$start$size, collapse = ";"),
               t(ci$checks), chain_mode_true_prob = chain_mode_prob[ch],
               chain_log_score = chain_scores[ch], stringsAsFactors = FALSE)
  }))

  list(summary = summary_row, chains = chain_table,
       diagnostics = cbind(setting = s_name, dataset = d, model = model_name, diag,
                           stringsAsFactors = FALSE),
       modes = data.frame(setting = s_name, dataset = d, model = model_name, location = locs,
                          prob_mode = p_all, prob_mode_given_true_count = p_cond,
                          stringsAsFactors = FALSE),
       scores = scores, chain_scores = chain_scores, tv = tv, tail_error = tail_error,
       p_true = mean(L_ind),
             priors = if (dist == "shifted_poisson") spois_priors else list(),
         y_train = y_train,
             plot = list(label = if (dist == "shifted_poisson") {
               sprintf("%s\nGamma(shape = %g, rate = %g)", model_name,
                 spois_priors$l0, spois_priors$L0)
             } else model_name,
             n_train = length(y_train),
             grid = grid, truth_p = truth_p, pbar = pbar,
                   band = band, y_test = y_test))
}

plot_curves <- function(path, s_name, d, fits) {
  panels <- length(fits)
  grDevices::png(path, width = max(1000L, 700L * panels), height = 800, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, panels), mar = c(4.4, 4.6, 3.6, 1), las = 1)
  reach <- function(grid, p) {
    i <- which(cumsum(p) >= 0.9995)
    if (length(i)) grid[i[1]] else max(grid)
  }
  for (p in fits) {
    grid <- p$grid
    xmax <- max(p$y_test, reach(grid, p$truth_p), reach(grid, p$pbar))
    visible <- grid <= xmax
    observed <- tabulate(as.integer(p$y_test) + 1L, nbins = length(grid)) / length(p$y_test)
    ymax <- 1.15 * max(observed[visible], p$band[2, visible], p$truth_p[visible])
    sample_label <- if (is.null(p$n_train)) "" else sprintf(" (n = %d)", p$n_train)
    graphics::plot(grid, p$truth_p, type = "n", xlim = c(0, xmax), ylim = c(0, ymax),
                   xlab = "Count", ylab = "Probability mass",
             main = sprintf("%s dataset %d%s: %s", s_name, d, sample_label, p$label))
    graphics::segments(grid, 0, grid, observed, col = "#CDD1D4", lwd = 4)
    graphics::polygon(c(grid, rev(grid)), c(p$band[1, ], rev(p$band[2, ])),
                      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.18), border = NA)
    graphics::lines(grid, p$truth_p, col = "#202020", lwd = 2)
    graphics::lines(grid, p$pbar, col = "#0072B2", lwd = 2, lty = 2)
  }
  graphics::legend("topright", c("True PMF", "Posterior mean PMF", "95% pointwise band", "Test frequencies"),
                   col = c("#202020", "#0072B2", "#BDD9E9", "#CDD1D4"),
                   lty = c(1, 2, 1, 1), lwd = c(2, 2, 8, 4), bty = "n", cex = 0.8)
}

results_path <- file.path(out_dir, "comparison-results.rds")
results <- if (reuse_other_models) {
  if (!file.exists(results_path)) stop("A selective rerun requires saved comparison results.")
  readRDS(results_path)
} else list()
refitted <- character()
for (s_index in setting_indices) {
  s_name <- names(settings)[s_index]
  s <- settings[[s_index]]
  for (d in 1:2) {
    y_train <- simulate_counts(s, n_train, 30000L + 10L * s_index + d)
    y_test <- simulate_counts(s, n_test, 31000L + 10L * s_index + d)
    if (reuse_other_models) {
      for (cached_model in setdiff(names(models), selected_model)) {
        cached_key <- paste(s_name, d, cached_model, sep = ".")
        cached <- results[[cached_key]]
        if (is.null(cached) || !identical(cached$plot$y_test, y_test)) {
          stop("Missing or incompatible saved comparison fit: ", cached_key)
        }
      }
    }
    for (m_index in seq_along(models)) {
      model_name <- names(models)[m_index]
      if (nzchar(selected_model) && model_name != selected_model) next
      dist <- models[[m_index]]
      key <- paste(s_name, d, model_name, sep = ".")
      results[[key]] <- NULL
      chains <- list()
      for (ch in seq_len(n_chains)) {
        if (proc.time()[["elapsed"]] - run_start > time_limit) break
        seed <- 50000L + 1000L * s_index + 100L * d + 10L * m_index + ch
        chains[[ch]] <- fit_chain(y_train, dist, seed)
        message(sprintf("%s chain %d done (%.1f s)", key, ch, chains[[ch]]$elapsed))
      }
      if (length(chains) < n_chains) {
        message(key, ": time limit reached after ", length(chains), " chains; not analysed")
        next
      }
      results[[key]] <- analyse_fit(s_name, s, d, model_name, dist, y_train, y_test, chains)
      refitted <- c(refitted, key)
      message(key, " analysed")
    }
  }
}

saveRDS(results, results_path)

for (s_name in names(settings)[setting_indices]) for (d in 1:2) {
  fits <- lapply(names(models), function(m) results[[paste(s_name, d, m, sep = ".")]]$plot)
  fits <- Filter(Negate(is.null), fits)
  if (length(fits)) {
    plot_curves(file.path(out_dir, sprintf("curves-%s-%d.png", s_name, d)), s_name, d, fits)
  }
}

paired <- list()
for (s_name in names(settings)) for (d in 1:2) {
  a <- results[[paste(s_name, d, "SNB", sep = ".")]]
  for (b_name in c("NB", "SPois")) {
    b <- results[[paste(s_name, d, b_name, sep = ".")]]
    if (is.null(a) || is.null(b)) next
    diff <- a$scores - b$scores
    paired[[length(paired) + 1]] <- data.frame(
      setting = s_name, dataset = d, comparison = paste("SNB minus", b_name),
      log_score_diff = mean(diff), paired_se_test_obs = sd(diff) / sqrt(length(diff)),
      chain_log_score_sd_SNB = sd(a$chain_scores), chain_log_score_sd_other = sd(b$chain_scores),
      tv_diff = a$tv - b$tv, tail_error_diff = a$tail_error - b$tail_error,
      prob_true_mode_count_diff = a$p_true - b$p_true, n_nonfinite_diff = sum(!is.finite(diff)),
      stringsAsFactors = FALSE)
  }
}

summary_table <- dplyr::bind_rows(lapply(results, `[[`, "summary"))
diagnostics <- dplyr::bind_rows(lapply(results, `[[`, "diagnostics"))
modes <- do.call(rbind, lapply(results, `[[`, "modes"))
chains_table <- dplyr::bind_rows(lapply(results, `[[`, "chains"))
paired_table <- do.call(rbind, paired)

utils::write.csv(summary_table, file.path(out_dir, "comparison-summary.csv"), row.names = FALSE)
if (length(paired)) {
  utils::write.csv(paired_table, file.path(out_dir, "comparison-paired.csv"), row.names = FALSE)
}
utils::write.csv(diagnostics, file.path(out_dir, "comparison-diagnostics.csv"), row.names = FALSE)
utils::write.csv(modes, file.path(out_dir, "comparison-modes.csv"), row.names = FALSE)
utils::write.csv(chains_table, file.path(out_dir, "comparison-chains.csv"), row.names = FALSE)

cat("\nTotal elapsed (s):", round(proc.time()[["elapsed"]] - run_start, 1), "\n")
cat("Fits sampled this run:", length(refitted), "\n")
expected_fits <- if (reuse_other_models) length(settings) * 2 * length(models) else {
  length(setting_indices) * 2 * if (nzchar(selected_model)) 1L else length(models)
}
cat("Fits analysed:", length(results), "of", expected_fits, "\n")
cat("Chain checks passed:", sum(chains_table$checks_pass), "of", nrow(chains_table), "\n")
cat("\nDiagnostic status counts:\n"); print(table(diagnostics$status))
cat("\nSummary:\n")
print(summary_table[, c("setting", "dataset", "model", "tv", "mean_log_score", "tail_error",
                        "prob_true_mode_count", "true_modes", "mean_pmf_modes", "mode_count_distribution")])
cat("\nPaired differences:\n")
print(paired_table)
