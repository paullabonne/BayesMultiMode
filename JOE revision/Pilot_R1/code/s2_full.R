# R1 comment 4a, full run. Usage: Rscript s2_full.R <workers>
# Shifted and unshifted NB, J = 10, 50,000 burn-in + 50,000 kept, every 25th draw summarized.
# Real data: 4 chains. Simulations: 50 data sets; chain 1 on every data set, 4 chains on the first 5.
# Unshifted NB with J = 20: G4 to G8, n = 2,000, 20 data sets.
# Finished tasks are skipped on a restart, so it can be restarted with more or fewer workers.
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
workers <- as.integer(commandArgs(TRUE)[1])
res10 <- file.path(S, "res_full"); res20 <- file.path(S, "res_full_J20")
dir.create(res10, showWarnings = FALSE); dir.create(res20, showWarnings = FALSE)

tasks <- list()
add <- function(design, n, rep, model, chain, J) tasks[[length(tasks) + 1]] <<- list(design = design, n = as.integer(n), rep = as.integer(rep), model = model, chain = as.integer(chain), J = J)
for (d in c("default", "dna", "inflation")) for (m in c("SNB", "NB")) for (ch in 1:4) add(d, 0, 1, m, ch, 10)
for (n in c(2000, 500)) for (r in 1:50) for (d in paste0("G", 1:8)) for (m in c("SNB", "NB"))
  for (ch in if (r <= 5) 1:4 else 1) add(d, n, r, m, ch, 10)
for (r in 1:20) for (d in paste0("G", 4:8)) add(d, 2000, r, "NB", 1, 20)
cat(length(tasks), "tasks,", workers, "workers, start", format(Sys.time()), "\n")

library(parallel)
cl <- makeCluster(workers)
clusterExport(cl, c("S", "res10", "res20"))
clusterEvalQ(cl, { source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); NULL })
out <- parLapplyLB(cl, tasks, function(tk) tryCatch(
  run_task(tk, if (tk$J == 10) res10 else res20, burnin = 50000, nkeep = 50000, J = tk$J, thin = 25),
  error = function(e) paste("ERROR", tk$design, tk$n, tk$rep, tk$model, tk$chain, conditionMessage(e))), chunk.size = 1)
stopCluster(cl)
err <- grep("^ERROR", unlist(out), value = TRUE)
cat("finished", format(Sys.time()), "| errors:", length(err), "\n"); if (length(err)) print(err)
