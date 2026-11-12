# R1 comment 4a: does an unshifted NB mixture fail where a shifted NB mixture succeeds?
# Designs, true probability mass functions, and the summary of one chain.

GRID <- 0:400

heap_y <- function(x, h) ifelse(runif(length(x)) < h, 5 * floor(x / 5 + 0.5), x)

# Each design: a simulator, the true pmf on GRID, and the true mode locations.
# G7 has tiny further modes above 20, each with probability below 0.003. They are not counted.
snb_pmf <- function(w, k, m, r) {
  p <- rep(0, length(GRID))
  for (j in seq_along(w)) {
    id <- GRID >= k[j]
    p[id] <- p[id] + w[j] * dnbinom(GRID[id] - k[j], size = r[j], mu = m[j])
  }
  p
}
heap_pmf <- function(m, r, h) {
  q <- dnbinom(0:2000, size = r, mu = m); p <- rep(0, length(GRID))
  for (i in seq_along(q)) {
    x <- i - 1; t <- 5 * floor(x / 5 + 0.5)
    if (x <= max(GRID)) p[x + 1] <- p[x + 1] + (1 - h) * q[i]
    if (t <= max(GRID)) p[t + 1] <- p[t + 1] + h * q[i]
  }
  p
}
rnorm_pmf <- function(w, m, s) {
  p <- rep(0, length(GRID))
  for (j in seq_along(w)) p <- p + w[j] * (pnorm(GRID + 0.5, m[j], s[j]) - pnorm(GRID - 0.5, m[j], s[j]))
  p / sum(p)
}

DESIGNS <- list(
  G1 = list(sim = function(n) simulate_snb(n, 1, 0, 5, 2),
            pmf = snb_pmf(1, 0, 5, 2), modes = 2, group = "control"),
  G2 = list(sim = function(n) ifelse(runif(n) < 0.6, 0, rnbinom(n, size = 4, mu = 7)),
            pmf = 0.6 * (GRID == 0) + 0.4 * dnbinom(GRID, size = 4, mu = 7), modes = c(0, 5), group = "control"),
  G3 = list(sim = function(n) simulate_snb(n, c(.6, .4), c(0, 0), c(3, 30), c(2, 5)),
            pmf = snb_pmf(c(.6, .4), c(0, 0), c(3, 30), c(2, 5)), modes = c(1, 23), group = "control"),
  G4 = list(sim = function(n) simulate_snb(n, c(.6, .4), c(0, 10), c(2, 5), c(5, 3)),
            pmf = snb_pmf(c(.6, .4), c(0, 10), c(2, 5), c(5, 3)), modes = c(1, 13), group = "shift"),
  G5 = list(sim = function(n) simulate_snb(n, c(.5, .5), c(0, 50), c(15, 1), c(10, 2)),
            pmf = snb_pmf(c(.5, .5), c(0, 50), c(15, 1), c(10, 2)), modes = c(13, 50), group = "shift"),
  G6 = list(sim = function(n) simulate_snb(n, c(.4, .3, .3), c(0, 7, 14), c(2, 3, 3), c(5, 5, 5)),
            pmf = snb_pmf(c(.4, .3, .3), c(0, 7, 14), c(2, 3, 3), c(5, 5, 5)), modes = c(1, 9, 16), group = "shift"),
  G7 = list(sim = function(n) heap_y(rnbinom(n, size = 8, mu = 8), 0.5),
            pmf = heap_pmf(8, 8, 0.5), modes = c(0, 5, 10, 15, 20), group = "other"),
  G8 = list(sim = function(n) {
              z <- runif(n) < 0.5
              pmax(0, round(ifelse(z, rnorm(n, 10, 1.5), rnorm(n, 25, 2))))
            },
            pmf = rnorm_pmf(c(.5, .5), c(10, 25), c(1.5, 2)), modes = c(10, 25), group = "other"))

# The limit for unshifted mixtures: for y >= 1, p(y) <= dpois(y, y), for any J.
unshifted_bound <- function(y) ifelse(y == 0, 1, dpois(y, y))
tv_floor <- function(p) sum(pmax(p - unshifted_bound(GRID), 0))

data_seed <- function(g, n, rep) 100000L + 1000L * g + 100L * (n == 2000) + rep

# Different starting values for each chain
chain_init <- function(y, J, chain) {
  if (chain == 1) return(list(rtilde = rep(1 / 11, J), kappa = rep(0, J)))   # r = 10
  top <- as.integer(quantile(y, 0.9))
  k <- if (top > 0) sample(0:top, J, replace = TRUE) else rep(0L, J)
  k[1] <- 0L                                                                  # keeps every y reachable
  list(rtilde = rep(1 / 3, J), kappa = k)                                     # r = 2
}

# Summary of one chain, from every 5th stored draw
summarize_chain <- function(o, y, true_modes, thin = 5) {
  idx <- seq(thin, nrow(o$mu), by = thin); K <- length(idx)
  tab <- table(y); u <- as.numeric(names(tab)); cnt <- as.numeric(tab)
  pm <- matrix(0, K, length(GRID)); L <- numeric(K); K0 <- numeric(K)
  pr_mode <- rep(0, length(GRID))
  at <- matrix(FALSE, K, length(true_modes)); near <- at
  for (s in seq_len(K)) {
    i <- idx[s]; act <- which(o$nj[i, ] > 0); w <- o$w[i, act] / sum(o$w[i, act])
    r <- (1 - o$rtilde[i, act]) / o$rtilde[i, act]
    for (a in seq_along(act)) {
      id <- GRID >= o$kappa[i, act[a]]
      pm[s, id] <- pm[s, id] + w[a] * dnbinom(GRID[id] - o$kappa[i, act[a]], size = r[a], mu = o$mu[i, act[a]])
    }
    m <- snb_modes(o$mu[i, ], o$rtilde[i, ], o$kappa[i, ], o$w[i, ], o$nj[i, ])
    L[s] <- m$L; K0[s] <- length(act)
    pr_mode[m$locs + 1] <- pr_mode[m$locs + 1] + 1
    for (t in seq_along(true_modes)) {
      at[s, t] <- true_modes[t] %in% m$locs
      near[s, t] <- any(abs(m$locs - true_modes[t]) <= 1)
    }
  }
  lpu <- log(pmax(pm[, u + 1, drop = FALSE], 1e-300))        # draws by unique values of y
  list(L = L, K0 = K0, loglik = as.numeric(lpu %*% cnt),
       pmf = colMeans(pm), pr_mode = pr_mode / K,
       at = colMeans(at), near = colMeans(near),
       lppd = sum(cnt * log(colMeans(exp(lpu)))),
       pwaic = sum(cnt * apply(lpu, 2, var)),
       acc = mean(o$acc[o$acc > 0]))
}

# Split Rhat on a set of chains, each a numeric vector of the same length
split_rhat <- function(chains) {
  h <- floor(length(chains[[1]]) / 2)
  s <- unlist(lapply(chains, function(x) list(x[1:h], x[(h + 1):(2 * h)])), recursive = FALSE)
  W <- mean(vapply(s, var, 0)); B <- h * var(vapply(s, mean, 0))
  sqrt(((h - 1) / h * W + B / h) / W)
}

# One task. Saved to its own file, so an interrupted run can resume.
real_data <- function(name) {
  if (name == "default") {
    x <- c(3002, 502, 187, 138, 233, 160, 107, 80, 59, 53, 41, 28, 34, 10, 13, 11, 4, 5, 8, 6, 3, 0, 1, 0,
           1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1)
    return(rep(0:34, x))
  }
  D <- Sys.getenv("PILOT_R1_DATA", unset = "C:/Users/JamieCross(LocalAdmi/Dropbox/Research/Papers/co_author/Herman/BayesMixturesShiftedNB/Code_corrected/MATLAB/Data/")   # data folder, with a trailing slash
  if (name == "dna") return(as.numeric(read.csv(paste0(D, "DATA_d4z4_4_411obs_a1_c2_y3.csv"))[[1]]))
  m <- read.csv(paste0(D, "michigan_15_april_2023.csv"))   # stand-in sample: March 2023, 0 to 50
  y <- m$PX1[m$YYYYMM == 202303]; y[!is.na(y) & y >= 0 & y <= 50]
}

run_task <- function(tk, res_dir, burnin = 5000, nkeep = 10000, J = 10, thin = 5) {
  f <- file.path(res_dir, sprintf("%s_n%d_r%d_%s_c%d.rds", tk$design, tk$n, tk$rep, tk$model, tk$chain))
  if (file.exists(f)) return(f)
  if (tk$design %in% names(DESIGNS)) {
    g <- match(tk$design, names(DESIGNS))
    set.seed(data_seed(g, tk$n, tk$rep)); y <- DESIGNS[[g]]$sim(tk$n); tm <- DESIGNS[[g]]$modes
  } else {
    y <- real_data(tk$design); tm <- integer(0)
  }
  set.seed(500000L + nchar(tk$design) * 7L + sum(utf8ToInt(tk$design)) * 11L + tk$n + 13L * tk$rep +
           (tk$model == "SNB") * 17L + 1009L * tk$chain)
  ini <- chain_init(y, J, tk$chain)
  t0 <- proc.time()[3]
  o <- snb_sampler(y, J, burnin, nkeep, sparse = TRUE, shifted = tk$model == "SNB", init = ini)
  secs <- proc.time()[3] - t0
  s <- summarize_chain(o, y, tm, thin = thin)
  s$secs <- secs; s$task <- tk; s$n_obs <- length(y)
  tmp <- paste0(f, ".tmp"); saveRDS(s, tmp); file.rename(tmp, f)
  f
}
