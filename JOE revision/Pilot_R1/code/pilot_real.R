# Pilot for the remaining R1 items on the three data sets: our model ("ours") and the shifted Poisson
# mixture ("sp") from NB_Comparisons.R, J = 10, sparse prior, 50,000 burn-in + 50,000 kept, 4 chains.
# Every kept draw is used for acceptance rates (r tilde, alpha); every 5th for modes, pmf, log-likelihood.
# Inflation: stand-in sample (March 2023, 0 to 50, n = 518); the paper's 480-observation sample is unknown.
# Usage: Rscript pilot_real.R <workers>
S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())   # working folder: set PILOT_R1_DIR or run from inside it
workers <- as.integer(commandArgs(TRUE)[1])
res <- file.path(S, "res_pilot"); dir.create(res, showWarnings = FALSE)

pilot_task <- function(tk) {
  f <- file.path(res, sprintf("%s_%s_c%d.rds", tk$data, tk$model, tk$chain)); if (file.exists(f)) return(f)
  y <- real_data(tk$data)
  seed <- 700000L + 101L * match(tk$data, c("default", "dna", "inflation")) + 1009L * tk$chain + (if (tk$model == "sp") 31L else 17L)
  set.seed(seed); ini <- chain_init(y, 10, tk$chain)
  o <- nb_sampler(y, J = 10, burnin = 50000, nkeep = 50000, algorithm = tk$model, thin = 1, init = ini, seed = seed)
  K <- nrow(o$mu)
  # Acceptance from consecutive kept draws: a non-empty component's rtilde changes only if accepted
  ne <- o$nj[-1, , drop = FALSE] > 0
  acc_rt <- if (tk$model == "sp") NA else sum((o$rtilde[-1, ] != o$rtilde[-K, ]) & ne) / sum(ne)
  acc_alpha <- mean(o$alpha[-1] != o$alpha[-K])
  # Every 5th draw: pmf on 0..max(y)+10, modes, log-likelihood
  idx <- seq(5, K, by = 5); G <- 0:(max(y) + 10)
  tab <- table(y); u <- as.numeric(names(tab)); cnt <- as.numeric(tab)
  pm <- matrix(0, length(idx), length(G)); L <- numeric(length(idx)); locs <- vector("list", length(idx))
  for (s in seq_along(idx)) {
    i <- idx[s]; act <- which(o$nj[i, ] > 0); w <- o$w[i, act] / sum(o$w[i, act])
    r <- (1 - o$rtilde[i, act]) / o$rtilde[i, act]
    top <- max(o$kappa[i, act] + pmax(0, floor(ifelse(is.infinite(r), o$mu[i, act], (r - 1) * o$mu[i, act] / r)))) + 1
    yg <- 0:max(top, max(G)); p <- rep(0, length(yg))
    for (k in seq_along(act)) { id <- yg >= o$kappa[i, act[k]]; p[id] <- p[id] + w[k] * dnbinom(yg[id] - o$kappa[i, act[k]], size = r[k], mu = o$mu[i, act[k]]) }
    pm[s, ] <- p[G + 1]
    q <- p[1:(top + 1)]; sg <- sign(diff(c(0, q))); nz <- which(sg != 0); lc <- integer(0)
    for (t in seq_along(nz)) if (sg[nz[t]] == 1 && (t == length(nz) || sg[nz[t + 1]] == -1)) lc <- c(lc, nz[t] - 1)
    L[s] <- length(lc); locs[[s]] <- lc
  }
  lpu <- log(pmax(pm[, u + 1, drop = FALSE], 1e-300)); ll <- as.numeric(lpu %*% cnt)
  act_w <- o$w * (o$nj > 0); act_w <- act_w / rowSums(act_w)
  disp <- if (tk$model == "sp") rep(NA, length(idx)) else rowSums(act_w[idx, ] * log((1 - o$rtilde[idx, ]) / o$rtilde[idx, ]))
  pr_mode <- rep(0, length(G)); for (lc in locs) pr_mode[lc[lc <= max(G)] + 1] <- pr_mode[lc[lc <= max(G)] + 1] + 1
  Ls <- sort(unique(L)); pr_mode_L <- sapply(Ls, function(l) { v <- rep(0, length(G)); id <- which(L == l)
    for (s in id) { lc <- locs[[s]]; v[lc[lc <= max(G)] + 1] <- v[lc[lc <= max(G)] + 1] + 1 }; v / length(id) })
  colnames(pr_mode_L) <- Ls
  out <- list(task = tk, n = length(y), G = G, cpu = o$cpu, acc_rt = acc_rt, acc_alpha = acc_alpha,
              sd_rt = o$sd1, alpha_step = 1, alpha = o$alpha[idx], loglik = ll, L = L, disp = disp,
              K0 = rowSums(o$nj[idx, ] > 0), pr_L = table(L) / length(L), pr_mode = pr_mode / length(idx), pr_mode_L = pr_mode_L,
              pmf = colMeans(pm), pmf_draws = pm[round(seq(1, nrow(pm), length.out = 200)), ],
              lppd = sum(cnt * log(colMeans(exp(lpu)))), pwaic = sum(cnt * apply(lpu, 2, var)))
  tmp <- paste0(f, ".tmp"); saveRDS(out, tmp); file.rename(tmp, f); f
}
tasks <- list()
for (d in c("default", "inflation", "dna")) for (m in c("ours", "sp")) for (ch in 1:4) tasks[[length(tasks) + 1]] <- list(data = d, model = m, chain = ch)
cat(length(tasks), "tasks,", workers, "workers, start", format(Sys.time()), "\n")
library(parallel)
cl <- makeCluster(workers)
clusterExport(cl, c("S", "res", "pilot_task"))
clusterEvalQ(cl, { source(file.path(S, "snb_s2_sampler.R")); source(file.path(S, "s2_functions.R")); source(file.path(S, "NB_Comparisons.R")); NULL })
out <- parLapplyLB(cl, tasks, function(tk) tryCatch(pilot_task(tk), error = function(e) paste("ERROR", tk$data, tk$model, tk$chain, conditionMessage(e))), chunk.size = 1)
stopCluster(cl)
err <- grep("^ERROR", unlist(out), value = TRUE)
cat("finished", format(Sys.time()), "| errors:", length(err), "\n"); if (length(err)) print(err)
