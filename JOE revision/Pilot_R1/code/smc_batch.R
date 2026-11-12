# Runs smc_ml.R for one model over the three data sets and three seeds, one after the other.
# Usage: Rscript smc_batch.R <model>
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
m <- commandArgs(TRUE)[1]; rs <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (seed in 1:3) for (d in c("dna", "inflation", "default")) {
  f <- file.path(S, "res_smc", sprintf("%s_%s_N1000_s5_seed%d.rds", d, m, seed))
  if (!file.exists(f)) system2(rs, c(shQuote(file.path(S, "smc_ml.R")), d, m, "1000", "5", seed))
}
cat("batch done", m, format(Sys.time()), "\n")
