phase2_settings <- function() {
  settings = sparse_prior_settings()
  settings$seed = 20261106L
  settings$replications = 100L
  settings$nb_iter = 12000L
  settings$burnin = 6000L
  settings$methods = c("neg_binomial", "shifted_poisson", "shifted_neg_binomial")
  settings$precision = list(proportion_mcse = 0.05, paired_tv_mcse = 0.005,
                            paired_log_score_mcse = 0.01)
  settings
}

phase2_scenarios <- function() {
  pilot = pilot_scenarios()
  list(
    poisson_control = pilot$poisson_control,
    overlapping_nb = pilot$overlapping_nb,
    shifted_poisson = list(family = "poisson", weight = c(0.5, 0.5), shift = c(0, 10),
                           mean = c(1.4, 3.4), size = c(NA_real_, NA_real_), expected_modes = c(1, 13)),
    shifted_nb = pilot$shifted_nb,
    weak_mode = list(family = "nb", weight = c(0.9, 0.1), shift = c(0, 6),
                     mean = c(1, 5), size = c(2, 2), expected_modes = c(0, 7)),
    structural_zero = list(family = "nb", weight = c(0.55, 0.45), shift = c(0, 0),
                           mean = c(0, 8), size = c(2, 3), expected_modes = c(0, 5))
  )
}

phase2_data <- function(settings, scenarios) {
  result = lapply(seq_along(scenarios), function(scenario_index) {
    lapply(seq_len(settings$replications), function(replication) {
      seed = settings$seed + scenario_index * 1000000L + replication * 1000L
      set.seed(seed)
      list(seed = seed, training = simulate_counts(scenarios[[scenario_index]], settings$sample_size),
           test = simulate_counts(scenarios[[scenario_index]], settings$test_size))
    })
  })
  names(result) = names(scenarios)
  result
}

paired_comparisons <- function(results) {
  rows = list()
  for (scenario in unique(results$scenario)) {
    proposed = results[results$scenario == scenario & results$method == "shifted_neg_binomial", ]
    for (method in c("neg_binomial", "shifted_poisson")) {
      comparator = results[results$scenario == scenario & results$method == method, ]
      paired = merge(proposed, comparator, by = "replication", suffixes = c("_proposed", "_comparator"))
      for (metric in intersect(c("tv", "log_predictive", "mode_brier"), names(results))) {
        difference = paired[[paste0(metric, "_proposed")]] - paired[[paste0(metric, "_comparator")]]
        finite = is.finite(difference)
        usable = length(difference) > 1 && all(finite)
        validated = usable && all(paired$all_diagnostics_pass_proposed %in% TRUE) &&
          all(paired$all_diagnostics_pass_comparator %in% TRUE)
        rows[[length(rows) + 1L]] = data.frame(
          scenario = scenario, comparator = method, metric = metric, pairs = nrow(paired),
          nonfinite_pairs = sum(!finite), diagnostic_failures = sum(
            !(paired$all_diagnostics_pass_proposed %in% TRUE & paired$all_diagnostics_pass_comparator %in% TRUE)),
          exploratory_mean_difference = if (usable) mean(difference) else NA_real_,
          exploratory_mcse = if (usable) sd(difference) / sqrt(length(difference)) else NA_real_,
          validated = validated)
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}