arguments = commandArgs(trailingOnly = TRUE)
stopifnot(all(arguments %in% c("--checks-only", "--smoke")))
script_argument = grep("^--file=", commandArgs(), value = TRUE)
study_directory = dirname(normalizePath(sub("^--file=", "", script_argument), winslash = "/"))
repository = normalizePath(file.path(study_directory, "../.."), winslash = "/")
source(file.path(study_directory, "phase1.R"))
settings = pilot_settings()
smoke = "--smoke" %in% arguments
if (smoke) {
  settings$replications = 1L
  settings$chains = 2L
  settings$nb_iter = 400L
  settings$burnin = 200L
}
output_directory = file.path(study_directory, "results", if (smoke) "smoke" else "phase1")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
gates = data.frame(gate = character(), status = character(), detail = character())
add_gate = function(gate, status, detail) {
  gates <<- rbind(gates, data.frame(gate, status, detail))
  write.csv(gates, file.path(output_directory, "gates.csv"), row.names = FALSE)
}
write_report = function() {
  report = c("# Phase 1 Status", "", paste("Run:", format(Sys.time(), tz = "UTC")),
             "", "Target: overall mixture PMF and its modes; no component recovery.", "",
             paste0("- ", gates$gate, ": ", gates$status, ". ", gates$detail), "",
             "Production is not authorised while any required gate is blocked, failed, or incomplete.",
             if (smoke) "This is a reduced smoke run, not the planned pilot." else "")
  writeLines(report, file.path(output_directory, "status.md"))
  writeLines(capture.output(sessionInfo()), file.path(output_directory, "session-info.txt"))
}
manifest = list(settings = settings,
                code_hashes = tools::md5sum(c(file.path(repository, "R/gibbs_sfm_algos.R"),
                                              file.path(repository, "tests/testthat/test-gibbs_sfm_algos.R"),
                                              file.path(study_directory, "phase1.R"),
                                              file.path(study_directory, "test-phase1.R"),
                                              file.path(study_directory, "run-phase1.R"))),
                commit = system2("git", c("-C", shQuote(repository), "rev-parse", "HEAD"), stdout = TRUE))
manifest_file = file.path(output_directory, "manifest.rds")
if (file.exists(manifest_file)) {
  previous = readRDS(manifest_file)
  if (!identical(previous, manifest)) stop("Existing results have different settings or code; archive them before rerunning.")
} else saveRDS(manifest, manifest_file)

checks = tryCatch({
  library(assertthat)
  source(file.path(repository, "R/gibbs_sfm_algos.R"))
  source(file.path(repository, "R/mixture_functions.R"))
  source(file.path(repository, "R/mix_mode.R"))
  testthat::test_file(file.path(repository, "tests/testthat/test-gibbs_sfm_algos.R"),
                      stop_on_failure = TRUE, reporter = "summary")
  testthat::test_file(file.path(study_directory, "test-phase1.R"),
                      stop_on_failure = TRUE, reporter = "summary")
  TRUE
}, error = function(error) {
  add_gate("correctness_checks", "failed", conditionMessage(error))
  FALSE
})
if (!checks) {
  write_report()
  stop("Phase 1 correctness checks failed.")
}
add_gate("correctness_checks", "passed", "Sampler regressions, independent PMF/mode truth, balance, priors, and initialisation checks passed.")

scenarios = pilot_scenarios()
truth = lapply(scenarios, scenario_truth, tolerance = settings$tail_tolerance)
saveRDS(truth, file.path(output_directory, "truth.rds"))
data_sets = lapply(seq_along(scenarios), function(scenario_index) {
  lapply(seq_len(settings$replications), function(replication) {
    seed = settings$seed + scenario_index * 100000L + replication * 1000L
    set.seed(seed)
    list(seed = seed, training = simulate_counts(scenarios[[scenario_index]], settings$sample_size),
         test = simulate_counts(scenarios[[scenario_index]], settings$test_size))
  })
})
names(data_sets) = names(scenarios)
saveRDS(data_sets, file.path(output_directory, "data.rds"))

observations = c(0, 1, 2, 3, 5, 8, 12)
reference = nb_quadrature(observations, settings$quadrature_size)
refined = nb_quadrature(observations, 2 * settings$quadrature_size)
baseline_file = file.path(output_directory, "one-component.rds")
if (file.exists(baseline_file)) {
  baseline = readRDS(baseline_file)
} else {
  baseline_settings = settings
  baseline_settings$nb_iter = 6000L
  baseline_settings$burnin = 2000L
  baseline = lapply(list(proposed = rneg_binom, traditional = traditional_nb_update), function(update) {
    lapply(seq_len(settings$chains), function(chain) {
      one_component_chain(observations, update, settings$seed + chain,
                           baseline_settings, if (identical(update, rneg_binom)) 0.1 else 1)
    })
  })
  saveRDS(baseline, baseline_file)
}
saveRDS(list(coarse = reference, refined = refined), file.path(output_directory, "quadrature.rds"))
baseline_check = lapply(baseline, function(chains) {
  estimates = do.call(rbind, lapply(chains, colMeans))
  list(mean = colMeans(estimates), chain_means = estimates,
       quadrature = refined, error = colMeans(estimates) - refined,
       acceptance = lapply(chains, attr, "acceptance"),
       proposal_scale = vapply(chains, attr, numeric(1), which = "proposal_scale"),
       adaptation_frozen = vapply(chains, attr, logical(1), which = "adaptation_frozen"))
})
saveRDS(baseline_check, file.path(output_directory, "baseline-check.rds"))
if (!all(vapply(baseline_check, function(check) all(check$adaptation_frozen), logical(1)))) {
  add_gate("equivalent_posterior_baseline", "failed", "Proposal adaptation continued after burn-in.")
  write_report()
  stop("Baseline adaptation gate failed.")
}

dependencies = trimws(strsplit(read.dcf(file.path(repository, "DESCRIPTION"))[1, "Imports"], ",")[[1]])
missing = c("pkgload", dependencies)[!vapply(c("pkgload", dependencies), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  add_gate("package_dependencies", "blocked", paste(missing, collapse = ", "))
  add_gate("equivalent_posterior_baseline", "incomplete", "Both coordinate samplers and independent quadrature saved; posterior diagnostics pending dependencies.")
  add_gate("pilot", "blocked", "No replacement sampler or regenerated data used.")
  write_report()
  stop("Missing package imports block the diagnostic and full pilot gates; see results status.md.")
}
pkgload::load_all(repository, quiet = TRUE)
testthat::test_local(repository, filter = "gibbs_sfm_algos|mix_mode", stop_on_failure = TRUE)
add_gate("package_dependencies", "passed", "All declared imports and pkgload available.")
baseline_summaries = lapply(baseline, summarise_chains)
baseline_pass = vapply(baseline_summaries, function(summary) {
  bounds = pmax(5 * summary$mcse_mean, 2 * abs(refined - reference), 0.002)
  diagnostics_pass(summary, settings) && all(abs(summary$mean - refined) <= bounds)
}, logical(1))
saveRDS(baseline_summaries, file.path(output_directory, "baseline-diagnostics.rds"))
write.csv(do.call(rbind, lapply(names(baseline_summaries), function(method) {
  cbind(method, baseline_summaries[[method]], quadrature = unname(refined))
})), file.path(output_directory, "baseline-diagnostics.csv"), row.names = FALSE)
if (!all(baseline_pass)) {
  add_gate("equivalent_posterior_baseline", "failed", "Coordinate samplers did not both meet quadrature agreement and R-hat/ESS gates.")
  add_gate("pilot", "blocked", "Baseline gate must pass first; no production run started.")
  write_report()
  stop("Equivalent-posterior baseline gate failed.")
}

traditional_fit = traditional_fit_function()
fixed_data = data_sets$overlapping_nb[[1]]$training
fixed_complexity = tryCatch(lapply(seq_len(settings$chains), function(chain) {
  set.seed(settings$seed + 500000L + chain)
  baseline_priors = pilot_priors("neg_binomial", settings)
  baseline_priors$sig = 1
  fit = traditional_fit(fixed_data, K = 2, dist = "neg_binomial", priors = baseline_priors,
                         nb_iter = settings$nb_iter, burnin = settings$burnin, print = FALSE)
  stopifnot(all(is.finite(fit$mcmc)))
  fit
}), error = function(error) error)
if (inherits(fixed_complexity, "error")) {
  add_gate("equivalent_posterior_baseline", "failed", conditionMessage(fixed_complexity))
  add_gate("pilot", "blocked", "Traditional mixture adapter feasibility must pass first.")
  write_report()
  stop("Traditional mixture feasibility gate failed.")
}
saveRDS(fixed_complexity, file.path(output_directory, "traditional-mixture.rds"))
add_gate("equivalent_posterior_baseline", "passed", "Both single-NB samplers agree with quadrature; traditional mixture adapter runs at fixed complexity. This is feasibility, not a speed claim.")
if ("--checks-only" %in% arguments) {
  add_gate("pilot", "incomplete", "Checks-only run; planned 480 pilot chains not started.")
  write_report()
  quit(status = 0)
}

methods = c("poisson", "neg_binomial", "shifted_poisson", "shifted_neg_binomial")
rows = list()
cores = parallel::detectCores(logical = FALSE)
workers = if (is.na(cores)) 1L else min(settings$workers, settings$chains, max(1L, cores - 1L))
cluster = if (workers > 1) parallel::makePSOCKcluster(workers, outfile = "") else NULL
tryCatch({
if (!is.null(cluster)) {
  parallel::clusterCall(cluster, function(repository, study_directory) {
    pkgload::load_all(repository, quiet = TRUE)
    source(file.path(study_directory, "phase1.R"))
    NULL
  }, repository, study_directory)
}
for (scenario_index in seq_along(scenarios)) {
  for (replication in seq_len(settings$replications)) {
    data_set = data_sets[[scenario_index]][[replication]]
    for (method_index in seq_along(methods)) {
      method = methods[method_index]
      key = paste(names(scenarios)[scenario_index], replication, method, sep = "-")
      result_file = file.path(output_directory, paste0(key, ".rds"))
      cat("Pilot:", key, "\n")
      if (file.exists(result_file)) {
        result = readRDS(result_file)
      } else {
        result = if (is.null(cluster)) {
          lapply(seq_len(settings$chains), pilot_chain, data_set = data_set,
                   method_index = method_index, method = method, settings = settings)
        } else {
          parallel::parLapply(cluster, seq_len(settings$chains), pilot_chain, data_set = data_set,
                              method_index = method_index, method = method, settings = settings)
        }
        saveRDS(result, result_file)
      }
      failures = vapply(result, function(chain) !is.na(chain$error), logical(1))
      row = data.frame(scenario = names(scenarios)[scenario_index], replication, method,
                        failures = sum(failures), passed = FALSE,
                        elapsed = sum(vapply(result, function(chain) chain$elapsed, numeric(1))))
      if (!any(failures)) {
        pooled = result[[1]]$fit
        pooled$mcmc = do.call(rbind, lapply(result, function(chain) chain$fit$mcmc))
        pooled$mcmc_all = pooled$mcmc
        pooled$loglik = unlist(lapply(result, function(chain) {
          tail(chain$fit$loglik, nrow(chain$fit$mcmc))
        }))
        evaluated = tryCatch(mixture_chain_summary(pooled, truth[[scenario_index]], data_set$test,
                         getFromNamespace("discrete_MF", "BayesMultiMode"), settings$tail_tolerance),
                         error = function(error) error)
        if (inherits(evaluated, "error")) {
          row$evaluation_error = conditionMessage(evaluated)
        } else {
          retained = nrow(result[[1]]$fit$mcmc)
          chains = lapply(seq_len(settings$chains), function(chain) {
            evaluated$diagnostics[(chain - 1) * retained + seq_len(retained), , drop = FALSE]
          })
          diagnostic = summarise_chains(chains)
          row$passed = diagnostics_pass(diagnostic, settings) &&
            all(is.finite(evaluated$metrics)) &&
            !any(lengths(lapply(result, function(chain) chain$warnings)))
          row = cbind(row, as.data.frame(as.list(c(evaluated$metrics,
                                    modal_metrics(evaluated$mode_locations, truth[[scenario_index]])))))
          saveRDS(list(diagnostics = diagnostic, metrics = row),
                   file.path(output_directory, paste0(key, "-summary.rds")))
        }
      }
      rows[[length(rows) + 1L]] = row
      columns = unique(unlist(lapply(rows, names)))
      table = do.call(rbind, lapply(rows, function(row) {
        row[setdiff(columns, names(row))] = NA
        row[columns]
      }))
      write.csv(table, file.path(output_directory, "pilot.csv"), row.names = FALSE)
    }
  }
}
}, finally = {
  if (!is.null(cluster)) parallel::stopCluster(cluster)
})
add_gate("pilot", if (smoke) "incomplete" else if (all(table$passed)) "passed" else "failed",
         paste(sum(table$passed), "of", nrow(table), "method/dataset groups passed diagnostic and warning gates."))
write_report()
if (!smoke && !all(table$passed)) stop("Pilot gates failed; inspect saved failures and diagnostics before production.")