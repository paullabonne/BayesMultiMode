arguments = commandArgs(trailingOnly = TRUE)
if (!length(arguments)) arguments = "--prepare"
stopifnot(length(arguments) == 1, arguments %in% c("--prepare", "--smoke", "--calibrate", "--run"))
mode = substring(arguments, 3)
script_argument = grep("^--file=", commandArgs(), value = TRUE)
study_directory = dirname(normalizePath(sub("^--file=", "", script_argument), winslash = "/"))
repository = normalizePath(file.path(study_directory, "../.."), winslash = "/")
source(file.path(study_directory, "phase1.R"))
source(file.path(study_directory, "phase2.R"))
settings = phase2_settings()
scenarios = phase2_scenarios()
truth = lapply(scenarios, scenario_truth, tolerance = settings$tail_tolerance)
files = c(list.files(file.path(repository, "R"), pattern = "[.]R$", full.names = TRUE),
          file.path(study_directory, c("phase1.R", "phase2.R", "run-phase2.R")))
protocol = list(settings = settings, scenarios = scenarios, code_hashes = tools::md5sum(files),
                commit = system2("git", c("-C", shQuote(repository), "rev-parse", "HEAD"), stdout = TRUE))
directory = file.path(study_directory, "results", "phase2", mode)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
manifest_file = file.path(directory, "manifest.rds")
if (file.exists(manifest_file)) {
  stopifnot(identical(readRDS(manifest_file), protocol))
} else saveRDS(protocol, manifest_file)
if (mode %in% c("smoke", "calibrate")) settings$replications = 1L
if (mode == "smoke") {
  settings$nb_iter = 400L
  settings$burnin = 200L
  settings$chains = 2L
}
if (mode == "run") {
  calibration_file = file.path(study_directory, "results", "phase2", "calibrate", "gate.rds")
  if (!file.exists(calibration_file)) stop("Run and inspect --calibrate before production.")
  gate = readRDS(calibration_file)
  if (!identical(gate$protocol, protocol) || !isTRUE(gate$passed)) {
    stop("Calibration has not passed under this protocol; production is not authorised.")
  }
}
data_sets = phase2_data(settings, scenarios)
saveRDS(data_sets, file.path(directory, "data.rds"))
saveRDS(truth, file.path(directory, "truth.rds"))
png(file.path(directory, "population-pmfs.png"), width = 1500, height = 1000, res = 150)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (scenario in names(truth)) {
  population = truth[[scenario]]
  upper = which(cumsum(population$pmf) >= 0.999)[1] - 1
  plot(population$support, population$pmf, type = "h", xlim = c(0, upper),
        xlab = "Count", ylab = "Probability", main = scenario)
  points(population$modes, population$pmf[population$modes + 1], col = "firebrick", pch = 19)
}
dev.off()
if (mode == "prepare") {
  writeLines(c("# Phase 2 Prepared", "", "Six scenarios and their population modes are independently checked.",
               "The protocol, 600 independent training/test datasets, and population PMF figure are saved.",
               "No models fitted. The next commands are --smoke and --calibrate; --run requires a passing calibration.",
               "The structural-zero model includes an exact point mass (NB with mean zero), outside the fitted interior parameter space."),
             file.path(directory, "status.md"))
  cat("Phase 2 prepared:", length(scenarios) * settings$replications, "datasets.\n")
  quit(status = 0)
}
pkgload::load_all(repository, quiet = TRUE)
cores = parallel::detectCores(logical = FALSE)
workers = if (is.na(cores)) 1L else min(settings$workers, settings$chains, max(1L, cores - 1L))
cluster = parallel::makePSOCKcluster(workers, outfile = "")
rows = list()
tryCatch({
  parallel::clusterCall(cluster, function(repository, study_directory) {
    pkgload::load_all(repository, quiet = TRUE)
    source(file.path(study_directory, "phase1.R"))
    NULL
  }, repository, study_directory)
  for (scenario in names(scenarios)) {
    for (replication in seq_len(settings$replications)) {
      for (method_index in seq_along(settings$methods)) {
        method = settings$methods[method_index]
        key = paste(scenario, replication, method, sep = "-")
        cat("Phase 2", mode, ":", key, "\n")
        result_file = file.path(directory, paste0(key, ".rds"))
        data_set = data_sets[[scenario]][[replication]]
        if (file.exists(result_file)) result = readRDS(result_file) else {
          result = parallel::parLapply(cluster, seq_len(settings$chains), pilot_chain,
                                      data_set = data_set, method_index = method_index,
                                      method = method, settings = settings)
          saveRDS(result, result_file)
        }
        summary = tryCatch(evaluate_chain_group(result, truth[[scenario]], data_set$test, settings),
                            error = function(error) list(row = data.frame(
                              failures = sum(vapply(result, function(chain) !is.na(chain$error), logical(1))),
                              all_diagnostics_pass = FALSE, evaluation_error = conditionMessage(error))))
        saveRDS(summary, file.path(directory, paste0(key, "-summary.rds")))
        rows[[length(rows) + 1L]] = cbind(data.frame(scenario, replication, method), summary$row)
        results = dplyr::bind_rows(rows)
        write.csv(results, file.path(directory, "results.csv"), row.names = FALSE)
      }
    }
  }
}, finally = parallel::stopCluster(cluster))
passed = mode != "smoke" && all(results$all_diagnostics_pass %in% TRUE)
saveRDS(list(protocol = protocol, passed = passed), file.path(directory, "gate.rds"))
write.csv(paired_comparisons(results), file.path(directory, "paired-comparisons.csv"), row.names = FALSE)
writeLines(c("# Phase 2 Status", "", paste("Mode:", mode),
             paste("Completed groups:", nrow(results), "; chains:", nrow(results) * settings$chains),
             paste("Fit errors:", sum(results$failures), "; diagnostic pass count:", sum(results$all_diagnostics_pass %in% TRUE)),
             "PMF and modal diagnostics are reported separately. Constant mode indicators have undefined MCSE, not zero uncertainty.",
             "No unsuccessful datasets were replaced or dropped. Paired summaries are explicitly exploratory unless all paired groups pass.",
             if (passed) "The diagnostic gate passed for this run." else "This run does not authorise production claims."),
           file.path(directory, "status.md"))
writeLines(capture.output(sessionInfo()), file.path(directory, "session-info.txt"))
if (mode != "smoke" && !passed) stop("Phase 2 diagnostic gate failed; inspect retained results and flags.")