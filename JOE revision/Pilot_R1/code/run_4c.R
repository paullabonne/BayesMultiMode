# Referee 1, comment 4c: mixture of shifted NB (ours) against mixture of shifted Poisson ("sp"), both from
# NB_Comparisons.R, on the eight designs of comment 4a plus two shifted Poisson designs (P1, P2).
# Settings as in 4a: J = 10, sparse prior, 50,000 burn-in + 50,000 kept, summaries every 25th draw;
# 50 data sets at n = 500 and n = 2,000 (chain 1), and chains 2-4 on data sets 1-5.
# Ours on G1-G8 is already in res_full (same algorithm, snb_s2_sampler.R); here it is run on P1 and P2.
# Usage: Rscript run_4c.R <workers>
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
workers <- as.integer(commandArgs(TRUE)[1])
res <- file.path(S, "res_4c"); dir.create(res, showWarnings = FALSE)

setup <- function() {
  source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R"))
  sp_pmf <- function(w, k, m) { p <- rep(0, length(GRID)); for (j in seq_along(w)) { id <- GRID >= k[j]; p[id] <- p[id] + w[j] * dpois(GRID[id] - k[j], m[j]) }; p }
  sp_sim <- function(n, w, k, m) { z <- sample.int(length(w), n, TRUE, prob = w); k[z] + rpois(n, m[z]) }
  DESIGNS$P1 <<- list(sim = function(n) sp_sim(n, c(.4, .3, .3), c(0, 7, 14), c(2.5, 2.5, 2.5)),
                      pmf = sp_pmf(c(.4, .3, .3), c(0, 7, 14), c(2.5, 2.5, 2.5)), modes = c(2, 9, 16), group = "poisson")
  DESIGNS$P2 <<- list(sim = function(n) sp_sim(n, c(.6, .4), c(0, 20), c(3.5, 4.5)),
                      pmf = sp_pmf(c(.6, .4), c(0, 20), c(3.5, 4.5)), modes = c(3, 24), group = "poisson")
  invisible(NULL)
}
run_4c_task <- function(tk) {
  f <- file.path(res, sprintf("%s_n%d_r%d_%s_c%d.rds", tk$design, tk$n, tk$rep, tk$model, tk$chain)); if (file.exists(f)) return(f)
  g <- match(tk$design, names(DESIGNS)); set.seed(data_seed(g, tk$n, tk$rep)); y <- DESIGNS[[g]]$sim(tk$n)
  seed <- 600000L + 7L * g + tk$n + 13L * tk$rep + 1009L * tk$chain + (if (tk$model == "sp") 31L else 17L)
  set.seed(seed); ini <- chain_init(y, 10, tk$chain)
  o <- nb_sampler(y, J = 10, burnin = 50000, nkeep = 50000, algorithm = tk$model, thin = 1, init = ini, seed = seed)
  o$acc <- o$acc1
  s <- summarize_chain(o, y, DESIGNS[[g]]$modes, thin = 25)
  s$cpu <- o$cpu; s$task <- tk; s$n_obs <- length(y)
  tp <- DESIGNS[[g]]$pmf; s$tv <- 0.5 * sum(abs(s$pmf - tp))
  s$elpd_true <- sum(tp[tp > 0] * log(pmax(s$pmf[tp > 0], 1e-300)))     # expected log score under the truth
  tmp <- paste0(f, ".tmp"); saveRDS(s, tmp); file.rename(tmp, f); f
}
setup()
tasks <- list()
for (n in c(2000L, 500L)) for (r in 1:50) for (d in c("P1", "P2", "G5", "G6", "G7", "G8", "G4", "G1", "G2", "G3")) {
  models <- if (d %in% c("P1", "P2")) c("sp", "ours") else "sp"
  chains <- if (r <= 5) 1:4 else 1L
  for (m in models) for (ch in chains) tasks[[length(tasks) + 1]] <- list(design = d, n = n, rep = r, model = m, chain = ch)
}
cat(length(tasks), "tasks,", workers, "workers, start", format(Sys.time()), "\n")
library(parallel)
cl <- makeCluster(workers)
clusterExport(cl, c("S", "res", "setup", "run_4c_task"))
clusterEvalQ(cl, { setup(); NULL })
out <- parLapplyLB(cl, tasks, function(tk) tryCatch(run_4c_task(tk), error = function(e) paste("ERROR", tk$design, tk$n, tk$rep, tk$model, tk$chain, conditionMessage(e))), chunk.size = 1)
stopCluster(cl)
err <- grep("^ERROR", unlist(out), value = TRUE)
cat("finished", format(Sys.time()), "| errors:", length(err), "\n"); if (length(err)) print(err)
