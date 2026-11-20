test_that("static Gaussian Galaxy fits reproduce the broad count distributions in Tables D.1 and D.2", {
  values <- BayesMultiMode::galaxy
  width <- diff(range(values))
  priors <- list(b0 = mean(range(values)), B0 = width^2, c0 = 2, g0 = 0.2, G0 = 10 / width^2)
  starts <- c(8L, 15L, 22L, 29L)
  chains <- lapply(seq_along(starts), function(chain) {
    set.seed(41000L + chain)
    fit <- bayes_fit(values, K = 30, dist = "normal", priors = priors, nb_iter = 1300, burnin = 300,
                     print = FALSE, sampler = "telescoping", K_init = starts[chain],
                     weight_model = "static", K_prior = "uniform", concentration = 1)
    fit$state[, c("K", "Kplus"), drop = FALSE]
  })
  draws <- do.call(rbind, chains)
  paper <- cbind(K = c(0, 0, .060, .135, .188, .195, .158, .109, .068, .039, .022, .012, .006, .003, .002),
                 Kplus = c(0, 0, .070, .161, .228, .228, .159, .087, .040, .017, .006, .002, .001, 0, 0))
  target <- apply(paper, 2, cumsum)
  fitted <- apply(draws, 2, function(count) cumsum(tabulate(count, nbins = 30))[1:15] / length(count))
  expect_lt(max(abs(fitted[, "K"] - target[, "K"])), .15)
  expect_lt(max(abs(fitted[, "Kplus"] - target[, "Kplus"])), .15)
  crossing <- mean(draws[, "Kplus"] <= 6 & draws[, "K"] > 6)
  expect_lt(abs(crossing - (target[6, "Kplus"] - target[6, "K"])), .06)
})