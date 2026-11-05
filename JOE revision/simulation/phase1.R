pilot_settings <- function() {
  list(seed = 20261105L, sample_size = 1000L, test_size = 2000L,
      replications = 10L, chains = 4L, workers = 4L, K = 5L,
       nb_iter = 3000L, burnin = 1500L, tail_tolerance = 1e-10,
       max_rhat = 1.01, min_ess = 400, quadrature_size = 160L,
       priors = list(a0 = 10, A0 = 50, l0 = 2, L0 = 1, sig = 0.1))
}

sparse_prior_settings <- function() {
  settings = pilot_settings()
  settings$priors$a0 = 1
  settings$priors$A0 = 200
  settings
}

pilot_priors <- function(distribution, settings) {
  fields = if (distribution %in% c("poisson", "shifted_poisson")) {
    c("a0", "A0", "l0", "L0")
  } else {
    stopifnot(distribution %in% c("neg_binomial", "shifted_neg_binomial"))
    c("a0", "A0", "sig")
  }
  settings$priors[fields]
}

pilot_chain <- function(chain, data_set, method_index, method, settings) {
  seed = data_set$seed + method_index * 100L + chain
  set.seed(seed)
  started = proc.time()["elapsed"]
  warnings = character()
  tryCatch({
    fit = withCallingHandlers(
      BayesMultiMode::bayes_fit(data_set$training, settings$K, method,
                                pilot_priors(method, settings), settings$nb_iter, settings$burnin, FALSE),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
    stopifnot(all(is.finite(fit$mcmc)), all(is.finite(fit$loglik)))
    list(seed = seed, fit = fit, elapsed = unname(proc.time()["elapsed"] - started),
         warnings = warnings, error = NA_character_)
  }, error = function(error) {
    list(seed = seed, elapsed = unname(proc.time()["elapsed"] - started),
         warnings = warnings, error = conditionMessage(error))
  })
}

pilot_scenarios <- function() {
  list(
    poisson_control = list(family = "poisson", weight = 1, shift = 0,
                           mean = 3.4, size = NA_real_, expected_modes = 3),
    overlapping_nb = list(family = "nb", weight = c(0.5, 0.5), shift = c(0, 0),
                          mean = c(4, 6), size = c(2, 2), expected_modes = 2),
    shifted_nb = list(family = "nb", weight = c(0.5, 0.5), shift = c(0, 10),
                      mean = c(1, 5), size = c(2, 2), expected_modes = c(0, 12))
  )
}

simulate_counts <- function(scenario, count) {
  allocation = sample.int(length(scenario$weight), count, replace = TRUE,
                          prob = scenario$weight)
  observations = numeric(count)
  for (component in seq_along(scenario$weight)) {
    selected = which(allocation == component)
    if (scenario$family == "poisson") {
      observations[selected] = scenario$shift[component] +
        rpois(length(selected), scenario$mean[component])
    } else {
      observations[selected] = scenario$shift[component] +
        rnbinom(length(selected), mu = scenario$mean[component],
                size = scenario$size[component])
    }
  }
  observations
}

log_sum_exp <- function(values) {
  maximum = max(values)
  if (maximum == -Inf) return(-Inf)
  if (!is.finite(maximum)) stop("Invalid log probabilities.")
  maximum + log(sum(exp(values - maximum)))
}

mixture_log_pmf <- function(observations, components) {
  probabilities = vapply(seq_along(components$weight), function(component) {
    shifted = observations - components$shift[component]
    if (components$family == "poisson") {
      density = dpois(shifted, components$mean[component], log = TRUE)
    } else {
      density = dnbinom(shifted, mu = components$mean[component],
                        size = components$size[component], log = TRUE)
    }
    log(components$weight[component]) + density
  }, numeric(length(observations)))
  probabilities = matrix(probabilities, nrow = length(observations))
  matrixStats::rowLogSumExps(probabilities)
}

certified_grid <- function(components, tolerance = 1e-10) {
  if (components$family == "poisson") {
    rightmost_modes = components$shift + floor(components$mean)
    tail_limits = components$shift + qpois(1 - tolerance, components$mean)
  } else {
    rightmost_modes = components$shift +
      floor(pmax(0, (components$size - 1) * components$mean / components$size))
    tail_limits = components$shift + qnbinom(1 - tolerance, mu = components$mean,
                                             size = components$size)
  }
  upper = max(rightmost_modes + 2, tail_limits)
  if (!is.finite(upper) || upper > 100000) stop("PMF grid exceeds the pilot budget.")
  0:as.integer(upper)
}

independent_modes <- function(probabilities, support, tolerance = 1e-12) {
  stopifnot(length(probabilities) == length(support),
            all(is.finite(probabilities)), all(probabilities >= 0),
            all(diff(support) == 1))
  differences = diff(probabilities)
  equal = abs(differences) <= tolerance *
    pmax(head(probabilities, -1), tail(probabilities, -1))
  starts = which(c(TRUE, !equal))
  ends = c(starts[-1] - 1L, length(probabilities))
  peaks = vapply(seq_along(starts), function(index) {
    first = starts[index]
    last = ends[index]
    left = if (first == 1) 0 else probabilities[first - 1]
    right = if (last == length(probabilities)) 0 else probabilities[last + 1]
    probabilities[first] > left * (1 + tolerance) &&
      probabilities[last] > right * (1 + tolerance)
  }, logical(1))
  support[starts[peaks]]
}

scenario_truth <- function(scenario, tolerance = 1e-10) {
  support = certified_grid(scenario, tolerance)
  probabilities = exp(mixture_log_pmf(support, scenario))
  modes = independent_modes(probabilities, support)
  stopifnot(identical(as.numeric(modes), as.numeric(scenario$expected_modes)),
            abs(sum(probabilities) - 1) <= tolerance * 2)
  list(support = support, pmf = probabilities, modes = modes, components = scenario,
       mean = sum(scenario$weight * (scenario$shift + scenario$mean)))
}

draw_components <- function(draw, distribution) {
  indices = as.integer(sub("^eta", "", grep("^eta[0-9]+$", names(draw), value = TRUE)))
  stopifnot(length(indices) > 0)
  poisson = distribution %in% c("poisson", "shifted_poisson")
  shifted = distribution %in% c("shifted_poisson", "shifted_neg_binomial")
  list(family = if (poisson) "poisson" else "nb",
       weight = unname(draw[paste0("eta", indices)]),
       shift = if (shifted) unname(draw[paste0("kappa", indices)]) else rep(0, length(indices)),
       mean = unname(draw[paste0(if (poisson) "lambda" else "mu", indices)]),
       size = if (poisson) rep(NA_real_, length(indices)) else
         1 / unname(draw[paste0("rtilde", indices)]) - 1)
}

package_mode_locations <- function(components, support, mode_function) {
  parameters = cbind(eta = components$weight, shift = components$shift,
                     mean = components$mean, size = components$size)
  density = function(observations, parameters) {
    if (components$family == "poisson") {
      dpois(observations - parameters["shift"], parameters["mean"])
    } else {
      dnbinom(observations - parameters["shift"], mu = parameters["mean"],
               size = parameters["size"])
    }
  }
  locations = withCallingHandlers(
    mode_function(parameters, density, range(support), type = "all"),
    warning = function(warning) {
      if (conditionMessage(warning) == "Some modes are flat.") invokeRestart("muffleWarning")
    })
  locations[c(TRUE, diff(locations) > 1)]
}

mean_predictive_log_pmf <- function(log_probabilities) {
  matrixStats::rowLogSumExps(log_probabilities) - log(ncol(log_probabilities))
}

traditional_nb_update <- function(Y, std, d_mu, d_rtilde, accept_ratio,
                                  burnin = FALSE, mu_max) {
  if (!length(Y)) return(c(runif(1, 0, mu_max), runif(1, 0, 0.5), 0, std))
  size = 1 / d_rtilde - 1
  latent = rgamma(length(Y), shape = Y + size, scale = d_mu / (size + d_mu))
  shape = size * length(Y) - 1
  scale = 1 / (size * sum(latent))
  log_tail = pgamma(1 / mu_max, shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  d_mu = 1 / qgamma(log_tail + log(runif(1)), shape, scale = scale,
                     lower.tail = FALSE, log.p = TRUE)
  if (burnin) {
    if (accept_ratio < 0.1) std = max(std * 0.8, 0.01)
    if (accept_ratio > 0.4) std = std * 1.2
  }
  candidate = rnorm(1, size, std)
  accepted = 0
  if (is.finite(candidate) && candidate > 1) {
    log_ratio = sum(dnbinom(Y, size = candidate, mu = d_mu, log = TRUE)) -
      sum(dnbinom(Y, size = size, mu = d_mu, log = TRUE)) -
      2 * log1p(candidate) + 2 * log1p(size)
    if (!is.finite(log_ratio)) stop("Invalid traditional-coordinate acceptance ratio.")
    if (log(runif(1)) < min(0, log_ratio)) {
      size = candidate
      accepted = 1
    }
  }
  c(d_mu, 1 / (1 + size), accepted, std)
}

one_component_chain <- function(observations, update, seed, settings, scale) {
  set.seed(seed)
  mean_value = runif(1, 0.2, 0.8) * max(observations)
  transformed_size = runif(1, 0.05, 0.45)
  accepted = 0
  retained_accepted = 0
  retained_scales = numeric(settings$nb_iter - settings$burnin)
  result = matrix(NA_real_, settings$nb_iter - settings$burnin, 3,
                   dimnames = list(NULL, c("mean", "pmf_zero", "tail_eight")))
  for (iteration in seq_len(settings$nb_iter)) {
    draw = update(observations, scale, mean_value, transformed_size,
                   accepted / iteration, iteration <= settings$burnin, max(observations))
    mean_value = draw[1]
    transformed_size = draw[2]
    accepted = accepted + draw[3]
    scale = draw[4]
    if (iteration > settings$burnin) {
      retained_accepted = retained_accepted + draw[3]
      retained_scales[iteration - settings$burnin] = scale
      size = 1 / transformed_size - 1
      result[iteration - settings$burnin, ] = c(
        mean_value, dnbinom(0, mu = mean_value, size = size),
        pnbinom(7, mu = mean_value, size = size, lower.tail = FALSE))
    }
  }
  attr(result, "acceptance") = c(burnin = (accepted - retained_accepted) / settings$burnin,
                                 retained = retained_accepted / nrow(result))
  attr(result, "proposal_scale") = scale
  attr(result, "adaptation_frozen") = length(unique(retained_scales)) == 1
  result
}

nb_quadrature <- function(observations, resolution) {
  mean_grid = (seq_len(resolution) - 0.5) * max(observations) / resolution
  transformed_grid = (seq_len(resolution) - 0.5) * 0.5 / resolution
  grid = expand.grid(mean = mean_grid, transformed_size = transformed_grid)
  size = 1 / grid$transformed_size - 1
  log_weight = vapply(seq_len(nrow(grid)), function(index) {
    sum(dnbinom(observations, mu = grid$mean[index], size = size[index], log = TRUE))
  }, numeric(1))
  weight = exp(log_weight - log_sum_exp(log_weight))
  c(mean = sum(weight * grid$mean),
    pmf_zero = sum(weight * dnbinom(0, mu = grid$mean, size = size)),
    tail_eight = sum(weight * pnbinom(7, mu = grid$mean, size = size, lower.tail = FALSE)))
}

mixture_chain_summary <- function(fit, truth, test_data, mode_function, tolerance) {
  draws = fit$mcmc
  components = lapply(seq_len(nrow(draws)), function(index) draw_components(draws[index, ], fit$dist))
  upper = max(truth$support, vapply(components, function(component) {
    max(certified_grid(component, tolerance))
  }, numeric(1)))
  support = 0:upper
  points = sort(unique(c(support, test_data)))
  log_pmf = vapply(components, function(component) mixture_log_pmf(points, component),
                    numeric(length(points)))
  log_predictive = mean_predictive_log_pmf(log_pmf)
  modes = lapply(components, package_mode_locations, support = support, mode_function = mode_function)
  mode_count = lengths(modes)
  truth_pmf = exp(mixture_log_pmf(support, truth$components))
  predicted_pmf = exp(log_predictive[match(support, points)])
  tail_error = max(0, 1 - sum(predicted_pmf)) + max(0, 1 - sum(truth_pmf))
  pmf_zero = exp(log_pmf[match(0, points), ])
  mixture_mean = vapply(components, function(component) {
    sum(component$weight * (component$shift + component$mean))
  }, numeric(1))
  mode_indices = vapply(modes, function(locations) {
    length(locations) == length(truth$modes) && all(locations == truth$modes)
  }, logical(1))
  retained = (nrow(fit$mcmc_all) - nrow(draws) + 1):nrow(fit$mcmc_all)
  stopifnot(length(fit$loglik[retained]) == nrow(draws))
  list(log_pmf = log_pmf, points = points, mode_locations = modes,
       diagnostics = cbind(mixture_mean, pmf_zero, mode_count,
                            true_mode_set = as.numeric(mode_indices)),
       metrics = c(tv = sum(abs(predicted_pmf - truth_pmf)) / 2,
                   tv_tail_bound = tail_error / 2,
                   log_predictive = mean(log_predictive[match(test_data, points)])))
}

traditional_fit_function <- function() {
  sampler = getFromNamespace("gibbs_SFM_neg_binomial", "BayesMultiMode")
  sampler_environment = new.env(parent = environment(sampler))
  sampler_environment$rneg_binom = traditional_nb_update
  environment(sampler) = sampler_environment
  fit = BayesMultiMode::bayes_fit
  fit_environment = new.env(parent = environment(fit))
  fit_environment$gibbs_SFM = function(y, K, nb_iter, nb_burnin, priors, print, dist) {
    stopifnot(dist == "neg_binomial")
    sampler(y, K, nb_iter, nb_burnin, priors, print)
  }
  environment(fit) = fit_environment
  fit
}

summarise_chains <- function(chains) {
  draws = simplify2array(chains)
  draws = aperm(draws, c(1, 3, 2))
  dimnames(draws)[[3]] = colnames(chains[[1]])
  summary = as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(draws), "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  summary$constant = vapply(seq_len(dim(draws)[3]), function(index) {
    length(unique(as.vector(draws[, , index]))) == 1
  }, logical(1))
  summary$discrete = vapply(seq_len(dim(draws)[3]), function(index) {
    all(draws[, , index] %% 1 == 0)
  }, logical(1))
  summary
}

diagnostics_pass <- function(summary, settings) {
  varying = summary[!summary$constant, , drop = FALSE]
  continuous = varying[!varying$discrete, , drop = FALSE]
  nrow(varying) > 0 && all(is.finite(varying$rhat)) &&
    all(varying$rhat <= settings$max_rhat) &&
    all(is.finite(varying$ess_bulk)) && all(varying$ess_bulk >= settings$min_ess) &&
    all(is.finite(continuous$ess_tail)) && all(continuous$ess_tail >= settings$min_ess)
}

evaluate_chain_group <- function(result, truth, test_data, settings) {
  failures = vapply(result, function(chain) !is.na(chain$error), logical(1))
  warnings = sum(lengths(lapply(result, function(chain) chain$warnings)))
  row = data.frame(failures = sum(failures), warnings = warnings,
                    pmf_diagnostics_pass = FALSE, modal_diagnostics_pass = FALSE,
                    all_diagnostics_pass = FALSE,
                    fit_seconds = sum(vapply(result, function(chain) chain$elapsed, numeric(1))))
  if (any(failures)) return(list(row = row, errors = lapply(result, function(chain) chain$error)))
  retained = vapply(result, function(chain) nrow(chain$fit$mcmc), integer(1))
  stopifnot(length(unique(retained)) == 1)
  pooled = result[[1]]$fit
  pooled$mcmc = do.call(rbind, lapply(result, function(chain) chain$fit$mcmc))
  pooled$mcmc_all = pooled$mcmc
  pooled$loglik = unlist(lapply(result, function(chain) tail(chain$fit$loglik, nrow(chain$fit$mcmc))))
  evaluated = mixture_chain_summary(pooled, truth, test_data,
                                    getFromNamespace("discrete_MF", "BayesMultiMode"), settings$tail_tolerance)
  chain_indices = lapply(seq_along(result), function(chain) {
    (chain - 1) * retained[chain] + seq_len(retained[chain])
  })
  chains = lapply(chain_indices, function(indices) evaluated$diagnostics[indices, , drop = FALSE])
  diagnostic = summarise_chains(chains)
  valid = warnings == 0 && all(is.finite(evaluated$metrics))
  row$pmf_diagnostics_pass = valid && diagnostics_pass(
    diagnostic[diagnostic$variable %in% c("mixture_mean", "pmf_zero"), ], settings)
  modal = diagnostic[diagnostic$variable %in% c("mode_count", "true_mode_set"), ]
  row$modal_diagnostics_pass = if (all(modal$constant)) NA else valid && diagnostics_pass(modal, settings)
  row$all_diagnostics_pass = valid && diagnostics_pass(diagnostic, settings)
  posterior_pmf = exp(evaluated$log_pmf)
  chain_pmf = vapply(chain_indices, function(indices) {
    rowMeans(posterior_pmf[, indices, drop = FALSE])
  }, numeric(nrow(posterior_pmf)))
  pair_tv = combn(seq_along(result), 2, function(pair) {
    sum(abs(chain_pmf[, pair[1]] - chain_pmf[, pair[2]])) / 2
  })
  row$max_between_chain_pmf_tv = max(pair_tv)
  mode_counts = lengths(evaluated$mode_locations)
  probabilities = lapply(0:max(mode_counts, length(truth$modes)), function(count) {
    indicator = matrix(as.numeric(mode_counts == count), ncol = length(result))
    constant = length(unique(as.vector(indicator))) == 1
    data.frame(mode_count = count, probability = mean(indicator),
               mcse = if (constant) NA_real_ else posterior::mcse_mean(indicator),
               constant = constant)
  })
  row = cbind(row, as.data.frame(as.list(c(evaluated$metrics,
                                          modal_metrics(evaluated$mode_locations, truth)))))
  list(row = row, diagnostics = diagnostic, mode_probabilities = do.call(rbind, probabilities),
       pmf = data.frame(count = evaluated$points, chain_pmf, check.names = FALSE))
}

modal_metrics <- function(mode_locations, truth) {
  counts = lengths(mode_locations)
  true_count = length(truth$modes)
  levels = 0:max(counts, true_count)
  probabilities = vapply(levels, function(count) mean(counts == count), numeric(1))
  ordering = order(-probabilities, levels)
  included = ordering[seq_len(which(cumsum(probabilities[ordering]) >= 0.95)[1])]
  correct = which(counts == true_count)
  location_error = if (length(correct)) mean(vapply(correct, function(index) {
    mean(abs(mode_locations[[index]] - truth$modes))
  }, numeric(1))) else NA_real_
  c(mode_brier = sum((probabilities - as.numeric(levels == true_count))^2),
    mode_coverage_95 = as.numeric(true_count %in% levels[included]),
    mode_set_size_95 = length(included),
    probability_true_count = mean(counts == true_count),
    probability_multimodal = mean(counts > 1),
    location_mae_given_count = location_error,
    missed_modes = mean(pmax(0, true_count - counts)),
    extra_modes = mean(pmax(0, counts - true_count)))
}