test_that("samplers sum log mixture probabilities over observations", {
  samplers = list(gibbs_SFM_normal, gibbs_SFM_skew_n, gibbs_SFM_poisson,
                  gibbs_SFM_shifted_poisson, gibbs_SFM_neg_binomial,
                  gibbs_SFM_shifted_neg_binomial)
  observations = rep(0:5, 2)
  expected_loglik = sum(dpois(observations, lambda = 3, log = TRUE))
  priors = list(a0 = 1, A0 = 200, b0 = 2, B0 = 25, c0 = 2.5,
                g0 = 0.5, G0 = 1, C0 = 1, D_xi = 1, D_psi = 1,
                l0 = 1.1, L0 = 1, sig = 0.1)

  for (sampler in samplers) {
    set.seed(123)
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$rdirichlet = function(count, alpha) matrix(c(0.6, 0.4), nrow = 1)
    sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
    sampler_env$draw_kap = function(...) list(0, 0)
    sampler_env$rtnorm = function(count, ...) rep(0, count)
    sampler_env$rmvnorm = function(count, mean, ...) matrix(mean, nrow = count)
    fixed_density = function(observations, ...) dpois(observations, lambda = 3)
    sampler_env$dnorm = fixed_density
    sampler_env$dsn = fixed_density
    sampler_env$dpois = fixed_density
    sampler_env$dnbinom = fixed_density

    arguments = list(y = observations, K = 2, nb_iter = 3, nb_burnin = 0,
                     priors = priors, print = FALSE)
    arguments = arguments[names(arguments) %in% names(formals(sampler))]
    result = do.call(sampler, arguments)

    expect_equal(unname(result[2:3, "loglik"]), rep(expected_loglik, 2))
  }
})

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
    result = update(observations, 0.1, 5, 0.2, 0.2, mu_max = 12)
    expect_equal(result[2:4], c(0.2, 0, 0.1))
    expect_equal(proposal_count, 1L)
  }

  candidate = 0.005
  result = update(observations, 0.1, 5, 0.2, 0.2, mu_max = 12)
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
  result = update(observations, 0.1, 5, 0.2, 0.2, mu_max = 5)
  expect_equal(result[2:3], c(0.4, 1))

  update_env$runif = function(...) (1 + acceptance_probability) / 2
  result = update(observations, 0.1, 5, 0.2, 0.2, mu_max = 5)
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
    result = update(observations, 0.1, 5, 0.2, acceptance_rate, mu_max = 12)
    expect_equal(result[4], 0.1)
  }

  result = update(observations, 0.1, 5, 0.2, 0, burnin = TRUE, mu_max = 12)
  expect_equal(result[4], 0.08)
  result = update(observations, 0.1, 5, 0.2, 0.5, burnin = TRUE, mu_max = 12)
  expect_equal(result[4], 0.12)
})

test_that("negative binomial mean draws use the fixed full-data bound", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$rgamma = function(n, ...) rep(1, n)
  update_env$rnorm = function(...) -1
  update_env$runif = function(...) 0.5

  zeros = update(c(0, 0), 0.1, 2, 0.2, 0.2, mu_max = 10)
  small_counts = update(c(0, 1), 0.1, 2, 0.2, 0.2, mu_max = 10)
  log_probability = pgamma(1 / 10, shape = 7, scale = 1 / 8,
                           lower.tail = FALSE, log.p = TRUE) + log(0.5)
  expected_mean = 1 / qgamma(log_probability, shape = 7, scale = 1 / 8,
                             lower.tail = FALSE, log.p = TRUE)
  expect_equal(zeros[1], expected_mean)
  expect_equal(small_counts[1], expected_mean)
  expect_gt(expected_mean, 1)
  expect_lt(expected_mean, 10)

  update_env$rgamma = function(n, ...) rep(1e5, n)
  result = update(c(0, 0), 0.1, 2, 0.2, 0.2, mu_max = 10)
  expect_gt(result[1], 0)
  expect_lte(result[1], 10)
})

test_that("empty negative binomial components draw directly from their priors", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  uniform_bounds = list()
  update_env$runif = function(n, min, max) {
    uniform_bounds[[length(uniform_bounds) + 1]] <<- c(min, max)
    (min + max) / 2
  }
  update_env$rgamma = function(...) stop("An empty component has no latent observations.")
  update_env$rnorm = function(...) stop("An empty component needs no MH proposal.")

  result = update(numeric(0), 0.1, 2, 0.2, 0, burnin = TRUE, mu_max = 10)
  expect_equal(result, c(5, 0.25, 0, 0.1))
  expect_equal(uniform_bounds, list(c(0, 10), c(0, 0.5)))
})

test_that("both negative binomial samplers pass fixed bounds and empty allocations", {
  for (sampler in list(gibbs_SFM_neg_binomial, gibbs_SFM_shifted_neg_binomial)) {
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$kmeans = function(y, ...) {
      list(cluster = rep(1L, length(y)), centers = matrix(c(0, 2), ncol = 1))
    }
    sampler_env$rdirichlet = function(n, alpha) matrix(alpha / sum(alpha), nrow = n)
    sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
    sampler_env$draw_kap = function(...) list(0, 0)
    sampler_env$sample.int = function(n, size) n
    component_counts = numeric(0)
    prior_bounds = numeric(0)
    initial_means = numeric(0)
    sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, accept_ratio, burnin, mu_max) {
      component_counts <<- c(component_counts, length(Y))
      prior_bounds <<- c(prior_bounds, mu_max)
      initial_means <<- c(initial_means, d_mu)
      rneg_binom(Y, std, d_mu, d_rtilde, accept_ratio, burnin, mu_max = mu_max)
    }

    result = sampler(c(0, 1, 2, 3, 4, 5), K = 2, nb_iter = 1, nb_burnin = 0,
                     priors = list(a0 = 1, A0 = 200, sig = 0.1), print = FALSE)
    expect_equal(component_counts, c(6, 0))
    expect_equal(prior_bounds, c(5, 5))
    expect_equal(initial_means, c(2.5, 2))
    expect_true(all(is.finite(result)))
    if ("kappa2" %in% colnames(result)) {
      expect_equal(unname(result[1, "kappa2"]), 5)
    }
  }
})