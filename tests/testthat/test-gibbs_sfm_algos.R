test_that("negative binomial Gaussian proposals are rejected outside support", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$qgamma = function(...) 1 / 5
  update_env$runif = function(...) 1e-12
  proposal_count = 0L
  candidate = 0.2
  update_env$rnorm = function(...) {
    proposal_count <<- proposal_count + 1L
    candidate
  }

  observations = c(0, 1, 2, 3, 5, 8, 12)
  for (candidate in c(-0.1, 0, 0.5, 0.6)) {
    proposal_count = 0L
    result = update(observations, 0.1, 5, 0.2, 0.2)
    expect_equal(result[2:4], c(0.2, 0, 0.1))
    expect_equal(proposal_count, 1L)
  }

  candidate = 0.005
  result = update(observations, 0.1, 5, 0.2, 0.2)
  expect_equal(result[2:3], c(0.005, 1))
})

test_that("negative binomial Gaussian acceptance follows the likelihood ratio", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$qgamma = function(...) 1 / 5
  update_env$rnorm = function(...) 0.4

  observations = rep(5, 10)
  current_loglik = sum(dnbinom(observations, size = 4, mu = 5, log = TRUE))
  candidate_loglik = sum(dnbinom(observations, size = 1.5, mu = 5, log = TRUE))
  acceptance_probability = exp(candidate_loglik - current_loglik)
  expect_lt(acceptance_probability, 1)

  update_env$runif = function(...) acceptance_probability / 2
  result = update(observations, 0.1, 5, 0.2, 0.2)
  expect_equal(result[2:3], c(0.4, 1))

  update_env$runif = function(...) (1 + acceptance_probability) / 2
  result = update(observations, 0.1, 5, 0.2, 0.2)
  expect_equal(result[2:3], c(0.2, 0))
})

test_that("negative binomial Gaussian proposal tuning stops after burn-in", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$qgamma = function(...) 1 / 5
  update_env$rnorm = function(...) 0.2
  observations = c(0, 1, 2, 3, 5, 8, 12)

  for (acceptance_rate in c(0, 0.5)) {
    result = update(observations, 0.1, 5, 0.2, acceptance_rate)
    expect_equal(result[4], 0.1)
  }

  result = update(observations, 0.1, 5, 0.2, 0, burnin = TRUE)
  expect_equal(result[4], 0.08)
  result = update(observations, 0.1, 5, 0.2, 0.5, burnin = TRUE)
  expect_equal(result[4], 0.12)
})