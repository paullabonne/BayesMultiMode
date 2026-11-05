script_argument = grep("^--file=", commandArgs(), value = TRUE)
arguments = commandArgs(trailingOnly = TRUE)
stopifnot(all(arguments %in% "--sparsity"))
sparsity_check = "--sparsity" %in% arguments
study_directory = dirname(normalizePath(sub("^--file=", "", script_argument), winslash = "/"))
repository = normalizePath(file.path(study_directory, "../.."), winslash = "/")
source(file.path(study_directory, "phase1.R"))
pkgload::load_all(repository, quiet = TRUE)
directory = file.path(study_directory, "results", "phase1")
output_directory = file.path(study_directory, "results", if (sparsity_check) "phase1-sparsity" else "phase1-diagnostics")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
manifest = readRDS(file.path(directory, "manifest.rds"))
settings = manifest$settings
if (!sparsity_check) {
  settings$nb_iter = 12000L
  settings$burnin = 6000L
}
probe_manifest = list(source_manifest = manifest, settings = settings,
                      script_hash = tools::md5sum(c(sub("^--file=", "", script_argument),
                                                    file.path(study_directory, "phase1.R"),
                                                    file.path(repository, "R/gibbs_sfm_algos.R"))),
                      sparsity_check = sparsity_check)
manifest_file = file.path(output_directory, "manifest.rds")
if (file.exists(manifest_file)) {
  stopifnot(identical(readRDS(manifest_file), probe_manifest))
} else saveRDS(probe_manifest, manifest_file)
data_sets = readRDS(file.path(directory, "data.rds"))
truth = readRDS(file.path(directory, "truth.rds"))
methods = c("poisson", "neg_binomial", "shifted_poisson", "shifted_neg_binomial")
if (sparsity_check) {
  comparison = expand.grid(scenario = names(data_sets), method = methods,
                            prior = c("pilot", "sparse"), stringsAsFactors = FALSE)
  cluster = parallel::makePSOCKcluster(min(settings$workers, settings$chains), outfile = "")
  rows = list()
  tryCatch({
    parallel::clusterCall(cluster, function(repository, study_directory) {
      pkgload::load_all(repository, quiet = TRUE)
      source(file.path(study_directory, "phase1.R"))
      NULL
    }, repository, study_directory)
    for (index in seq_len(nrow(comparison))) {
      probe = comparison[index, ]
      current_settings = settings
      if (probe$prior == "sparse") {
        current_settings$priors$a0 = 1
        current_settings$priors$A0 = 200
      }
      data_set = data_sets[[probe$scenario]][[1]]
      key = paste(probe$scenario, probe$method, probe$prior, sep = "-")
      cat("Sparsity check:", key, "\n")
      result_file = file.path(output_directory, paste0(key, ".rds"))
      if (file.exists(result_file)) result = readRDS(result_file) else {
        result = parallel::parLapply(cluster, seq_len(settings$chains), pilot_chain,
                                    data_set = data_set, method_index = match(probe$method, methods),
                                    method = probe$method, settings = current_settings)
        saveRDS(result, result_file)
      }
      summary = tryCatch(evaluate_chain_group(result, truth[[probe$scenario]], data_set$test, current_settings),
                          error = function(error) list(row = data.frame(evaluation_error = conditionMessage(error))))
      saveRDS(summary, file.path(output_directory, paste0(key, "-summary.rds")))
      rows[[index]] = cbind(probe, replication = 1L, summary$row)
      columns = unique(unlist(lapply(rows, names)))
      table = do.call(rbind, lapply(rows, function(row) {
        row[setdiff(columns, names(row))] = NA
        row[columns]
      }))
      write.csv(table, file.path(output_directory, "sparsity-comparison.csv"), row.names = FALSE)
    }
  }, finally = parallel::stopCluster(cluster))
  writeLines(capture.output(sessionInfo()), file.path(output_directory, "session-info.txt"))
  print(table[, intersect(c("scenario", "method", "prior", "failures", "warnings", "pmf_diagnostics_pass",
                             "modal_diagnostics_pass", "all_diagnostics_pass", "max_between_chain_pmf_tv",
                             "evaluation_error"), names(table))])
  quit(status = 0)
}
probes = data.frame(scenario = c("shifted_nb", "shifted_nb"), replication = c(4L, 7L),
                    method = c("shifted_poisson", "shifted_neg_binomial"))
cluster = parallel::makePSOCKcluster(min(settings$workers, settings$chains), outfile = "")
tryCatch({
parallel::clusterCall(cluster, function(repository, study_directory) {
  pkgload::load_all(repository, quiet = TRUE)
  source(file.path(study_directory, "phase1.R"))
  NULL
}, repository, study_directory)
rows = list()
for (index in seq_len(nrow(probes))) {
  probe = probes[index, ]
  key = paste(probe$scenario, probe$replication, probe$method, sep = "-")
  data_set = data_sets[[probe$scenario]][[probe$replication]]
  result_file = file.path(output_directory, paste0(key, ".rds"))
  cat("Long-chain check:", key, "\n")
  if (file.exists(result_file)) {
    result = readRDS(result_file)
  } else {
    result = parallel::parLapply(cluster, seq_len(settings$chains), pilot_chain,
                                data_set = data_set, method_index = match(probe$method, methods),
                                method = probe$method, settings = settings)
    saveRDS(result, result_file)
  }
  stopifnot(all(vapply(result, function(chain) is.na(chain$error), logical(1))))
  pooled = result[[1]]$fit
  pooled$mcmc = do.call(rbind, lapply(result, function(chain) chain$fit$mcmc))
  pooled$mcmc_all = pooled$mcmc
  pooled$loglik = unlist(lapply(result, function(chain) tail(chain$fit$loglik, nrow(chain$fit$mcmc))))
  evaluated = mixture_chain_summary(pooled, truth[[probe$scenario]], data_set$test,
                                    getFromNamespace("discrete_MF", "BayesMultiMode"), settings$tail_tolerance)
  retained = nrow(result[[1]]$fit$mcmc)
  chains = lapply(seq_len(settings$chains), function(chain) {
    evaluated$diagnostics[(chain - 1) * retained + seq_len(retained), , drop = FALSE]
  })
  diagnostic = summarise_chains(chains)
  original = readRDS(file.path(directory, paste0(key, "-summary.rds")))$diagnostics
  row = data.frame(group = key, iterations = settings$nb_iter,
                    passed = diagnostics_pass(diagnostic, settings) && all(is.finite(evaluated$metrics)),
                    max_rhat = max(diagnostic$rhat, na.rm = TRUE),
                    min_bulk_ess = min(diagnostic$ess_bulk, na.rm = TRUE),
                    original_max_rhat = max(original$rhat, na.rm = TRUE),
                    original_min_bulk_ess = min(original$ess_bulk, na.rm = TRUE))
  rows[[index]] = row
  saveRDS(list(diagnostics = diagnostic, original_diagnostics = original,
                metrics = evaluated$metrics, settings = settings),
           file.path(output_directory, paste0(key, "-summary.rds")))
  pdf(file.path(output_directory, paste0(key, "-traces.pdf")), width = 10, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (variable in seq_len(ncol(chains[[1]]))) {
    matplot(do.call(cbind, lapply(chains, function(chain) chain[, variable])), type = "l",
             lty = 1, col = seq_along(chains), xlab = "Retained iteration",
             ylab = colnames(chains[[1]])[variable], main = key)
  }
  dev.off()
  write.csv(do.call(rbind, rows), file.path(output_directory, "long-chain-checks.csv"), row.names = FALSE)
  print(row)
}
}, finally = parallel::stopCluster(cluster))

pilot = read.csv(file.path(directory, "pilot.csv"))
files = list.files(directory, pattern = "-summary[.]rds$", full.names = TRUE)
diagnostics = do.call(rbind, lapply(files, function(file) {
  cbind(group = sub("-summary.rds$", "", basename(file)), readRDS(file)$diagnostics)
}))
diagnostics$failed = !diagnostics$constant &
  (!is.finite(diagnostics$rhat) | diagnostics$rhat > settings$max_rhat |
   !is.finite(diagnostics$ess_bulk) | diagnostics$ess_bulk < settings$min_ess |
   (!diagnostics$discrete & (!is.finite(diagnostics$ess_tail) | diagnostics$ess_tail < settings$min_ess)))
failures = aggregate(failed ~ variable, diagnostics, sum)
write.csv(failures, file.path(output_directory, "diagnostic-failures.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_directory, "all-diagnostics.csv"), row.names = FALSE)
report = c("# Phase 1 Pilot Assessment", "",
           "Target: overall mixture PMF and its modes, not individual components.", "",
           paste(nrow(pilot) * manifest$settings$chains, "chains completed across", nrow(pilot), "method/dataset groups."),
           paste("Fitting errors:", sum(pilot$failures), "; groups passing the predeclared gates:", sum(pilot$passed), "."),
           "Correctness, package dependencies, and equivalent-posterior baseline gates passed.", "",
           "## Convergence", "",
           paste0("- ", failures$variable, ": ", failures$failed, " diagnostic failures."), "",
           "## Targeted Longer Chains", "",
           "These post-pilot checks reuse two selected failing groups, with unchanged data, seeds, priors, and thresholds.",
           "They use 12000 iterations and 6000 burn-in; they do not replace or silently discard any original pilot group.", "",
           paste0("- ", vapply(rows, function(row) paste(row$group, ": gate", if (row$passed) "passed" else "failed",
                  "; maximum R-hat", round(row$original_max_rhat, 3), "->", round(row$max_rhat, 3),
                  "; minimum bulk ESS", round(row$original_min_bulk_ess), "->", round(row$min_bulk_ess)), character(1))), "",
           "## Decision", "",
           "The original 3000-iteration pilot is executed, but not convergence-validated. Production remains blocked.",
           "Do not use its aggregate performance estimates as evidence of superiority: many groups have low modal/zero-mass ESS.",
           "Use the trace plots and longer-chain checks to set the next local tuning or mixing investigation; do not rerun the full grid blindly.")
writeLines(report, file.path(directory, "assessment.md"))
writeLines(capture.output(sessionInfo()), file.path(output_directory, "session-info.txt"))