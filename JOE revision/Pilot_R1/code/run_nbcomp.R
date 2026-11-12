# All five algorithms of NB_Comparisons.R on the stage 1 data: G1, G4, G5, G6 (n = 2,000, data sets 1-5,
# chain 1) and the default data (chains 1-4). Within a task, the algorithms run one after the other in the
# same R process, on the same data, with the same starting values and the same seed.
# Usage: Rscript run_nbcomp.R <workers>
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
workers <- as.integer(commandArgs(TRUE)[1])
res <- file.path(S, "res_nbcomp"); dir.create(res, showWarnings = FALSE)
ALG <- c("ours", "mur_rw", "pr_rw", "zhou", "zhou_mh")
run_nbcomp_task <- function(tk) {
  f <- file.path(res, sprintf("%s_r%d_c%d.rds", tk$design, tk$rep, tk$chain)); if (file.exists(f)) return(f)
  if (tk$design %in% names(DESIGNS)) { g <- match(tk$design, names(DESIGNS)); set.seed(data_seed(g, 2000, tk$rep)); y <- DESIGNS[[g]]$sim(2000); n <- 2000L }
  else { y <- real_data(tk$design); n <- 0L }
  seed <- 500000L + nchar(tk$design) * 7L + sum(utf8ToInt(tk$design)) * 11L + n + 13L * tk$rep + 17L + 1009L * tk$chain
  set.seed(seed); ini <- chain_init(y, 10, tk$chain)
  out <- lapply(setNames(ALG, ALG), function(a) {
    o <- nb_sampler(y, J = 10, burnin = 5000, nkeep = 10000, algorithm = a, thin = 5, init = ini, seed = seed)
    s <- nb_summary(o, y); s$task <- tk; s
  })
  tmp <- paste0(f, ".tmp"); saveRDS(out, tmp); file.rename(tmp, f); f
}
tasks <- list()
for (ch in 1:4) tasks[[length(tasks) + 1]] <- list(design = "default", rep = 1L, chain = ch)
for (r in 1:5) for (d in c("G5", "G6", "G4", "G1")) tasks[[length(tasks) + 1]] <- list(design = d, rep = r, chain = 1L)
cat(length(tasks), "tasks x", length(ALG), "algorithms,", workers, "workers, start", format(Sys.time()), "\n")
library(parallel)
cl <- makeCluster(workers)
clusterExport(cl, c("S", "res", "ALG", "run_nbcomp_task"))
# snb_s2_sampler.R supplies simulate_snb for the designs; NB_Comparisons.R is sourced last, so its functions are used
clusterEvalQ(cl, { source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R")); NULL })
out <- parLapplyLB(cl, tasks, function(tk) tryCatch(run_nbcomp_task(tk), error = function(e) paste("ERROR", tk$design, tk$rep, tk$chain, conditionMessage(e))), chunk.size = 1)
stopCluster(cl)
err <- grep("^ERROR", unlist(out), value = TRUE)
cat("finished", format(Sys.time()), "| errors:", length(err), "\n"); if (length(err)) print(err)
