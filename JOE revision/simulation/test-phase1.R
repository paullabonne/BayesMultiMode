testthat::test_that("the sparse-prior follow-up changes only the concentration prior", {
  original = pilot_settings()
  sparse = sparse_prior_settings()
  testthat::expect_equal(sparse$priors[c("a0", "A0")], list(a0 = 1, A0 = 200))
  sparse$priors[c("a0", "A0")] = original$priors[c("a0", "A0")]
  testthat::expect_identical(sparse, original)
})

testthat::test_that("each method receives only its supported pilot priors", {
  settings = pilot_settings()
  for (distribution in c("poisson", "shifted_poisson", "neg_binomial", "shifted_neg_binomial")) {
    priors = pilot_priors(distribution, settings)
    testthat::expect_warning(check_priors(priors, distribution, c(0, 1, 2, 3)), NA)
    testthat::expect_equal(priors[c("a0", "A0")], settings$priors[c("a0", "A0")])
  }
})

testthat::test_that("group evaluation separates PMF stability and constant modal indicators", {
  testthat::skip_if_not_installed("BayesMultiMode")
  settings = pilot_settings()
  settings$nb_iter = 60L
  settings$burnin = 30L
  settings$sample_size = 100L
  set.seed(51)
  data_set = list(seed = 51L, training = simulate_counts(pilot_scenarios()$poisson_control, 100))
  fits = lapply(1:2, pilot_chain, data_set = data_set, method_index = 1,
                 method = "poisson", settings = settings)
  summary = evaluate_chain_group(fits, scenario_truth(pilot_scenarios()$poisson_control), 0:15, settings)
  testthat::expect_equal(summary$row$failures, 0)
  testthat::expect_true(is.finite(summary$row$max_between_chain_pmf_tv))
  testthat::expect_equal(sum(summary$mode_probabilities$probability), 1)
  testthat::expect_true(all(is.na(summary$mode_probabilities$mcse[summary$mode_probabilities$constant])))
})

testthat::test_that("pilot truth is known on certified integer grids", {
  for (scenario in pilot_scenarios()) {
    truth = scenario_truth(scenario)
    testthat::expect_equal(truth$modes, scenario$expected_modes)
    testthat::expect_lt(abs(sum(truth$pmf) - 1), 2e-10)
  }
})

testthat::test_that("truth detection uses leftmost plateaus and zero boundaries", {
  testthat::expect_equal(independent_modes(c(4, 3, 1, 2, 3, 2) / 15, 0:5), c(0, 4))
  testthat::expect_equal(independent_modes(c(4, 4, 1, 3, 3, 2) / 17, 0:5), c(0, 3))
  testthat::expect_equal(independent_modes(c(0.6, 0.3, 0.1), 0:2), 0)
})

testthat::test_that("predictive scores average probabilities rather than log scores", {
  probabilities = matrix(c(0.2, 0.8, 0.8, 0.2), nrow = 2)
  testthat::expect_equal(mean_predictive_log_pmf(log(probabilities)), rep(log(0.5), 2))
  testthat::expect_equal(mean_predictive_log_pmf(matrix(-Inf, 1, 2)), -Inf)
  testthat::expect_equal(log_sum_exp(c(-1000, -1001)), -1000 + log1p(exp(-1)))
})

testthat::test_that("relabelled or split components represent the same target", {
  scenario = pilot_scenarios()$shifted_nb
  support = certified_grid(scenario)
  permuted = scenario
  for (field in c("weight", "shift", "mean", "size")) permuted[[field]] = rev(scenario[[field]])
  testthat::expect_equal(mixture_log_pmf(support, scenario), mixture_log_pmf(support, permuted))
  split = scenario
  split$weight = c(0.25, 0.25, 0.5)
  split$shift = c(0, 0, 10)
  split$mean = c(1, 1, 5)
  split$size = c(2, 2, 2)
  testthat::expect_equal(mixture_log_pmf(support, scenario), mixture_log_pmf(support, split))
})

testthat::test_that("public count generators reproduce random mixture moments", {
  scenario = pilot_scenarios()$shifted_nb
  set.seed(37)
  observations = simulate_counts(scenario, 50000)
  testthat::expect_true(all(observations >= 0 & observations %% 1 == 0))
  testthat::expect_equal(mean(observations), scenario_truth(scenario)$mean, tolerance = 0.02)
  set.seed(37)
  testthat::expect_equal(observations, simulate_counts(scenario, 50000))
})

testthat::test_that("independent truth agrees with package discrete mode detection", {
  for (scenario in pilot_scenarios()) {
    truth = scenario_truth(scenario)
    testthat::expect_equal(package_mode_locations(scenario, truth$support, discrete_MF), truth$modes)
  }
})

testthat::test_that("shift kernels satisfy finite-grid detailed balance", {
  observations = c(3, 4, 6)
  for (distribution in c("shifted_poisson", "shifted_neg_binomial")) {
    shift_grid = 0:3
    extra = if (distribution == "shifted_poisson") 5 else c(0.2, 5)
    log_target = vapply(shift_grid, function(shift) {
      if (distribution == "shifted_poisson") {
        sum(dpois(observations - shift, 5, log = TRUE))
      } else {
        sum(dnbinom(observations - shift, size = 4, mu = 5, log = TRUE))
      }
    }, numeric(1))
    target = exp(log_target - log_sum_exp(log_target))
    transition = matrix(0, 4, 4)
    update = draw_kap
    update_environment = new.env(parent = environment(update))
    environment(update) = update_environment
    step = 0
    uniform = 0
    update_environment$sample = function(...) step
    update_environment$runif = function(...) uniform
    for (current in shift_grid) {
      for (step in c(-1, 1)) {
        candidate = current + step
        if (candidate < 0 || candidate > 3) next
        acceptance = min(1, exp(log_target[candidate + 1] - log_target[current + 1]))
        uniform = acceptance / 2
        testthat::expect_equal(update(observations, c(current, extra), 0, 3, distribution)[[1]], candidate)
        if (acceptance < 1) {
          uniform = (1 + acceptance) / 2
          testthat::expect_equal(update(observations, c(current, extra), 0, 3, distribution)[[1]], current)
        }
        transition[current + 1, candidate + 1] = acceptance / 3
      }
      transition[current + 1, current + 1] = 1 - sum(transition[current + 1, ])
    }
    flow = target * transition
    testthat::expect_equal(flow, t(flow), tolerance = 1e-12)
  }
})

testthat::test_that("concentration log target matches an independent reference", {
  weights = c(0.7, 0.2, 0.1)
  for (concentration in c(0.005, 0.2, 2)) {
    expected = dgamma(concentration, 2, scale = 0.5, log = TRUE) +
      lgamma(3 * concentration) - 3 * lgamma(concentration) +
      sum((concentration - 1) * log(weights))
    testthat::expect_equal(post_e0(concentration, 2, 0.5, weights, log = TRUE), expected)
  }
  testthat::expect_error(post_e0(0.2, 2, 0.5, c(1, 0), log = TRUE), "strictly positive")
})

testthat::test_that("traditional coordinates use the transformed size prior", {
  update = traditional_nb_update
  update_environment = new.env(parent = environment(update))
  environment(update) = update_environment
  update_environment$rgamma = function(count, ...) rep(1, count)
  update_environment$qgamma = function(...) 1 / 5
  update_environment$rnorm = function(...) 8
  observations = c(0, 1, 2, 3, 5, 8)
  probability = min(1, exp(sum(dnbinom(observations, 8, mu = 5, log = TRUE)) -
                            sum(dnbinom(observations, 4, mu = 5, log = TRUE)) -
                            2 * log(9) + 2 * log(5)))
  update_environment$runif = function(...) probability / 2
  testthat::expect_equal(update(observations, 1, 5, 0.2, 0.2, mu_max = 8)[2:3], c(1 / 9, 1))
  update_environment$rnorm = function(...) 0.5
  testthat::expect_equal(update(observations, 1, 5, 0.2, 0.2, mu_max = 8)[2:3], c(0.2, 0))
})

testthat::test_that("count and location metrics never match components", {
  truth = list(modes = c(0, 12))
  metrics = modal_metrics(list(c(0, 12), c(0, 13), 0), truth)
  testthat::expect_equal(unname(metrics["probability_true_count"]), 2 / 3)
  testthat::expect_equal(unname(metrics["location_mae_given_count"]), 0.25)
  testthat::expect_equal(unname(metrics["missed_modes"]), 1 / 3)
  testthat::expect_equal(unname(metrics["mode_brier"]), 2 / 9)
})

testthat::test_that("posterior quadrature is stable under refinement", {
  observations = c(0, 1, 2, 3, 5, 8, 12)
  coarse = nb_quadrature(observations, 40)
  fine = nb_quadrature(observations, 80)
  testthat::expect_equal(coarse, fine, tolerance = 0.002)
})

testthat::test_that("unsupported initialisation fails rather than changing K", {
  for (sampler in list(gibbs_SFM_poisson, gibbs_SFM_shifted_poisson,
                       gibbs_SFM_neg_binomial, gibbs_SFM_shifted_neg_binomial)) {
    arguments = list(y = rep(0:1, 5), K = 3, nb_iter = 2, nb_burnin = 1,
                     priors = list(a0 = 10, A0 = 50, l0 = 2, L0 = 1, sig = 0.1), print = FALSE)
    arguments = arguments[names(arguments) %in% names(formals(sampler))]
    testthat::expect_error(do.call(sampler, arguments), "distinct")
  }
  testthat::expect_error(gibbs_SFM_neg_binomial(rep(0, 5), 2, 2, 1,
                       list(a0 = 10, A0 = 50, sig = 0.1), FALSE), "positive full-data maximum")
})

testthat::test_that("discrete diagnostics require bulk ESS without undefined tail ESS", {
  summary = data.frame(constant = c(FALSE, FALSE, TRUE), discrete = c(FALSE, TRUE, TRUE),
                       rhat = c(1, 1, NA), ess_bulk = c(1000, 1000, NA),
                       ess_tail = c(1000, NA, NA))
  testthat::expect_true(diagnostics_pass(summary, pilot_settings()))
  summary$ess_bulk[2] = 20
  testthat::expect_false(diagnostics_pass(summary, pilot_settings()))
  summary$ess_bulk[2] = 1000
  summary$ess_tail[1] = NA
  testthat::expect_false(diagnostics_pass(summary, pilot_settings()))
})

testthat::test_that("Phase 2 truth and seeds are fixed independently of fitted outcomes", {
  source("phase2.R", local = TRUE)
  scenarios = phase2_scenarios()
  testthat::expect_length(scenarios, 6)
  for (scenario in scenarios) testthat::expect_equal(scenario_truth(scenario)$modes, scenario$expected_modes)
  settings = phase2_settings()
  testthat::expect_equal(settings$priors[c("a0", "A0")], list(a0 = 1, A0 = 200))
  settings$replications = 2L
  first = phase2_data(settings, scenarios)
  testthat::expect_identical(first, phase2_data(settings, scenarios))
  testthat::expect_false(identical(first[[1]][[1]]$training, first[[1]][[2]]$training))
  zero = scenarios$structural_zero
  testthat::expect_equal(exp(mixture_log_pmf(0, zero)), 0.55 + 0.45 * dnbinom(0, mu = 8, size = 3))
})

testthat::test_that("paired comparisons retain failures instead of selecting successful fits", {
  source("phase2.R", local = TRUE)
  results = expand.grid(scenario = "case", replication = 1:2,
                        method = c("neg_binomial", "shifted_poisson", "shifted_neg_binomial"))
  results$all_diagnostics_pass = TRUE
  results$tv = c(0.1, 0.2, 0.1, 0.2, 0.05, 0.1)
  comparison = paired_comparisons(results)
  testthat::expect_true(all(comparison$validated))
  results$all_diagnostics_pass[1] = FALSE
  comparison = paired_comparisons(results)
  testthat::expect_false(comparison$validated[comparison$comparator == "neg_binomial"])
  testthat::expect_true(all(comparison$pairs == 2))
  results$tv[1] = NA_real_
  comparison = paired_comparisons(results)
  testthat::expect_true(is.na(comparison$exploratory_mean_difference[comparison$comparator == "neg_binomial"]))
})