test_that("single-component count samplers preserve allocation matrix dimensions", {
  for (distribution in c("poisson", "shifted_poisson", "neg_binomial", "shifted_neg_binomial")) {
    set.seed(12)
    priors = if (distribution %in% c("poisson", "shifted_poisson")) {
      list(a0 = 1, A0 = 200, l0 = 2, L0 = 1)
    } else list(a0 = 1, A0 = 200, sig = 0.1)
    draws = gibbs_SFM(c(2, 3, 4, 5, 7, 10, 14), 1, 20, 10, priors, FALSE, distribution)
    expect_equal(nrow(draws), 20)
    expect_true(all(is.finite(draws)))
    expect_equal(unname(draws[, "eta1"]), rep(1, 20))
  }
})

test_that("collapsed concentration target integrates Dirichlet weights", {
  for (concentration in c(0.005, 0.2, 2)) {
    counts = c(2, 1)
    integral = integrate(function(weight) {
      weight^2 * (1 - weight) * dbeta(weight, concentration, concentration)
    }, 0, 1, rel.tol = 1e-9)$value
    expected = dgamma(concentration, 1, scale = 1/200, log = TRUE) + log(integral)
    expect_equal(post_e0(concentration, 1, 1/200, log = TRUE, counts = counts),
                 expected, tolerance = 1e-7)
  }
  expect_true(is.finite(post_e0(0.005, 1, 1/200, log = TRUE, counts = c(1000, 0, 0, 0, 0))))
  expect_equal(post_e0(0.005, 1, 1/200, log = TRUE, counts = rep(0, 5)),
               dgamma(0.005, 1, scale = 1/200, log = TRUE))
})

test_that("all samplers update collapsed concentration before refreshing weights", {
  samplers = list(normal = gibbs_SFM_normal, skew_normal = gibbs_SFM_skew_n,
                   poisson = gibbs_SFM_poisson, shifted_poisson = gibbs_SFM_shifted_poisson,
                   neg_binomial = gibbs_SFM_neg_binomial,
                   shifted_neg_binomial = gibbs_SFM_shifted_neg_binomial)
  for (distribution in names(samplers)) {
    sampler = samplers[[distribution]]
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    updates = character()
    sampler_env$draw_e0 = function(e0, nu0, S0, p = NULL, counts = NULL) {
      expect_null(p)
      expect_equal(sum(counts), 6)
      updates <<- c(updates, "concentration")
      list(0.2, 1)
    }
    sampler_env$rdirichlet = function(n, alpha) {
      updates <<- c(updates, "weights")
      expect_equal(sum(alpha), 6.4)
      matrix(alpha / sum(alpha), nrow = n)
    }
    arguments = list(y = 0:5, K = 2, nb_iter = 3, nb_burnin = 1,
             priors = check_priors(list(), distribution, 0:5), print = FALSE)
    arguments = arguments[names(arguments) %in% names(formals(sampler))]
    set.seed(17)
    do.call(sampler, arguments)
    expect_equal(updates, rep(c("concentration", "weights"), 3 - (distribution == "normal")))
  }
})

test_that("continuous samplers tolerate underflowed weights for empty components", {
  observations = rep(seq(-2, 2, length.out = 20), 2)
  for (distribution in c("normal", "skew_normal")) {
    sampler = if (distribution == "normal") gibbs_SFM_normal else gibbs_SFM_skew_n
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$rdirichlet = function(n, alpha) {
      matrix(c(1, rep(0, length(alpha) - 1)), nrow = n)
    }
    set.seed(123)
    draws = sampler(observations, K = 10, nb_iter = 20,
                    priors = check_priors(list(), distribution, observations), print = FALSE)
    retained = draws[-1, , drop = FALSE]
    expect_true(all(is.finite(retained)))
    expect_true(all(retained[, "eta1"] == 1))
    expect_true(all(retained[, paste0("eta", 2:10)] == 0))
  }
})

test_that("shifted Poisson shift likelihood includes every observation", {
  observations = c(3, 4, 6, 7)
  parameters = c(2, 5)
  expect_equal(post_kap(observations, parameters, "shifted_poisson"),
               sum(dpois(observations - 2, lambda = 5, log = TRUE)))
  expect_equal(post_kap(rep(1000, 20), c(0, 1000), "shifted_poisson"),
               sum(dpois(rep(1000, 20), lambda = 1000, log = TRUE)))
})

test_that("concentration updates use a finite log ratio", {
  update = draw_e0
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  candidate = 0.3
  update_env$rnorm = function(...) candidate - 0.2
  update_env$runif = function(...) 0
  weights = c(0.9, rep(0.1 / 999, 999))
  result = update(0.2, 2, 0.5, weights)
  expect_equal(result, list(0.3, 1))

  candidate = -0.1
  expect_equal(update(0.2, 2, 0.5, weights), list(0.2, 0))
})

test_that("shift updates reject candidates outside conditional support", {
  update = draw_kap
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  step = -1
  update_env$sample = function(...) step
  update_env$runif = function(...) 0
  observations = c(2, 4, 5)
  expect_equal(update(observations, c(0, 3), 0, 2, "shifted_poisson"), list(0, 0))
  step = 1
  expect_equal(update(observations, c(2, 3), 0, 2, "shifted_poisson"), list(2, 0))
  expect_equal(update(observations, c(0, 3), 0, 2, "shifted_poisson"), list(1, 1))
  expect_error(update(observations, c(3, 3), 0, 2, "shifted_poisson"),
               "outside its conditional support")
})

test_that("empty shifted Poisson components refresh shifts from their prior", {
  sampler = gibbs_SFM_shifted_poisson
  sampler_env = new.env(parent = environment(sampler))
  environment(sampler) = sampler_env
  sampler_env$kmeans = function(y, ...) {
    list(cluster = rep(1L, length(y)), centers = matrix(c(1, 3), ncol = 1))
  }
  sampler_env$rdirichlet = function(n, alpha) matrix(alpha / sum(alpha), nrow = n)
  sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
  sampler_env$draw_kap = function(...) list(0, 0)
  prior_bound = NA_real_
  sampler_env$sample.int = function(n, size) {
    prior_bound <<- n
    n
  }
  result = sampler(0:5, 2, 1, list(a0 = 1, A0 = 200, l0 = 2, L0 = 1), FALSE)
  expect_equal(prior_bound, 6)
  expect_equal(unname(result[1, "kappa2"]), 5)
  expect_true(all(is.finite(result)))
})

test_that("joint shift-mean draws preserve the mean and use the conditional likelihood", {
  update = draw_kap_mu
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  observations = c(4, 4, 6, 7)
  shifts = 0:4
  means = 6 - shifts
  log_probabilities = vapply(seq_along(shifts), function(candidate) {
    sum(dnbinom(observations - shifts[candidate], size = 3,
                mu = means[candidate], log = TRUE))
  }, numeric(1))
  update_env$sample.int = function(n, size, prob) {
    expect_equal(n, length(shifts))
    expect_equal(size, 1)
    expect_equal(prob / sum(prob), exp(log_probabilities - max(log_probabilities)) /
                   sum(exp(log_probabilities - max(log_probabilities))))
    n
  }
  result = update(observations, kappa = 1, mu = 5, r_tilde = 0.25, mu_max = 7)
  expect_equal(result, list(kappa = 4, mu = 2))
  expect_equal(result$kappa + result$mu, 6)
})

test_that("joint shift-mean draws respect both shift and mean support", {
  update = draw_kap_mu
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$sample.int = function(n, size, prob) {
    expect_equal(n, 3)
    1L
  }
  expect_equal(update(c(8, 9), 7, 9, 0.25, mu_max = 10), list(kappa = 6, mu = 10))

  update_env$sample.int = function(n, size, prob) {
    expect_equal(n, 4)
    n
  }
  expect_equal(update(c(4, 5), 2, 2, 0.25, mu_max = 5), list(kappa = 3, mu = 1))

  update_env$sample.int = function(...) stop("A singleton support needs no random draw.")
  expect_equal(update(c(5, 6), 5, 6, 0.25, mu_max = 6), list(kappa = 5, mu = 6))
  expect_equal(update(c(0, 1), 0, 2, 0.25, mu_max = 5), list(kappa = 0, mu = 2))
  expect_error(update(c(0, 1), 0, 10, 0.25, mu_max = 5), "conditional support")
})

test_that("joint shift-mean probabilities remain stable for small likelihoods", {
  update = draw_kap_mu
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  observations = rep(c(4, 6, 7), 1000)
  expect_equal(exp(sum(dnbinom(observations - 1, size = 3, mu = 5, log = TRUE))), 0)
  update_env$sample.int = function(n, size, prob) {
    expect_true(all(is.finite(prob)))
    expect_equal(max(prob), 1)
    which.max(prob)
  }
  result = update(observations, 1, 5, 0.25, mu_max = 7)
  expect_equal(result$kappa + result$mu, 6)
  update_env$dnbinom = function(x, ...) rep(-Inf, length(x))
  expect_error(update(observations, 1, 5, 0.25, mu_max = 7), "Non-finite joint shift-mean")
})

test_that("shifted negative binomial sampler passes both joint draws to the next update", {
  sampler = gibbs_SFM_shifted_neg_binomial
  sampler_env = new.env(parent = environment(sampler))
  environment(sampler) = sampler_env
  sampler_env$kmeans = function(y, ...) {
    list(cluster = rep(1L, length(y)), centers = matrix(6, ncol = 1))
  }
  sampler_env$draw_kap = function(...) stop("The shift-only update must not be used.")
  sampler_env$draw_kap_mu = function(x, kappa, mu, r_tilde, mu_max) {
    expect_equal(x, c(5, 6, 8))
    expect_equal(c(kappa, mu, r_tilde, mu_max), c(0, 6, 0.4, 8))
    list(kappa = 4, mu = 2)
  }
  sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, mu_max) {
    expect_equal(Y, c(1, 2, 4))
    expect_equal(d_mu, 2)
    expect_equal(mu_max, 8)
    c(d_mu, d_rtilde, 0, std)
  }
  result = sampler(c(5, 6, 8), 1, 1, 0, list(a0 = 1, A0 = 200, sig = 0.1), FALSE)
  expect_equal(unname(result[1, c("kappa1", "mu1")]), c(4, 2))
})

test_that("joint shifted negative binomial fits retain valid draws and log likelihoods", {
  set.seed(117)
  observations = c(5 + rnbinom(80, size = 3, mu = 10),
                    50 + rnbinom(80, size = 2, mu = 1))
  fit = bayes_fit(observations, K = 2, dist = "shifted_neg_binomial",
                  nb_iter = 120, burnin = 75, print = FALSE)
  expect_equal(dim(fit$mcmc_all), c(120, 8))
  expect_equal(nrow(fit$mcmc), 45)
  expect_identical(colnames(fit$mcmc_all),
                   c("eta1", "eta2", "kappa1", "kappa2", "mu1", "mu2", "rtilde1", "rtilde2"))
  expect_true(all(is.finite(fit$mcmc_all)))
  weights = fit$mcmc_all[, c("eta1", "eta2")]
  shifts = fit$mcmc_all[, c("kappa1", "kappa2")]
  means = fit$mcmc_all[, c("mu1", "mu2")]
  shapes = fit$mcmc_all[, c("rtilde1", "rtilde2")]
  expect_equal(unname(rowSums(weights)), rep(1, 120))
  expect_true(all(weights >= 0))
  expect_true(all(shifts >= 0 & shifts <= max(observations) & shifts %% 1 == 0))
  expect_true(all(means > 0 & means <= max(observations)))
  expect_true(all(shapes > 0 & shapes < 0.5))
  expect_true(all(shifts[, 1] <= shifts[, 2]))
  for (draw in c(1, 50, 51, 75, 120)) {
    probabilities = vapply(1:2, function(component) {
      weights[draw, component] * dnbinom(observations - shifts[draw, component],
        size = 1 / shapes[draw, component] - 1, mu = means[draw, component])
    }, numeric(length(observations)))
    expect_equal(fit$loglik[draw], sum(log(rowSums(probabilities))), tolerance = 1e-8)
  }
})

test_that("count samplers never replace invalid allocations with uniform draws", {
  for (sampler in list(gibbs_SFM_poisson, gibbs_SFM_shifted_poisson,
                       gibbs_SFM_neg_binomial, gibbs_SFM_shifted_neg_binomial)) {
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$rdirichlet = function(n, alpha) matrix(c(0.5, 0.5), nrow = n)
    sampler_env$draw_kap = function(...) list(0, 0)
    sampler_env$draw_kap_mu = function(x, kappa, mu, ...) list(kappa = 0, mu = mu)
    sampler_env$rneg_binom = function(...) c(2, 0.2, 0, 0.1)
    zero_density = function(observations, ...) rep(0, length(observations))
    sampler_env$dpois = zero_density
    sampler_env$dnbinom = zero_density
    arguments = list(y = 0:5, K = 2, nb_iter = 2, nb_burnin = 1,
                     priors = list(a0 = 10, A0 = 50, l0 = 2, L0 = 1, sig = 0.1), print = FALSE)
    arguments = arguments[names(arguments) %in% names(formals(sampler))]
    expect_error(do.call(sampler, arguments), "Non-finite count allocation")
  }
})

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
    sampler_env$draw_kap_mu = function(x, kappa, mu, ...) list(kappa = 0, mu = mu)
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
    result = update(observations, 0.1, 5, 0.2, mu_max = 12)
    expect_equal(result[2:4], c(0.2, 0, 0.1))
    expect_equal(proposal_count, 1L)
  }

  candidate = 0.005
  result = update(observations, 0.1, 5, 0.2, mu_max = 12)
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
  result = update(observations, 0.1, 5, 0.2, mu_max = 5)
  expect_equal(result[2:3], c(0.4, 1))

  update_env$runif = function(...) (1 + acceptance_probability) / 2
  result = update(observations, 0.1, 5, 0.2, mu_max = 5)
  expect_equal(result[2:3], c(0.2, 0))
})

test_that("rtilde batch adaptation targets 30 to 40 percent and diminishes", {
  expect_equal(adapt_rtilde_proposals(rep(0.1, 6), c(0, 15, 17, 20, 50, 0),
                                     c(rep(50, 5), 0), 1),
               0.1 * exp(c(-0.1, 0, 0, 0, 0.1, 0)))
  expect_equal(adapt_rtilde_proposals(c(0.1, 0.1), c(0, 50), c(50, 50), 400),
               0.1 * exp(c(-0.05, 0.05)))
  expect_lt(adapt_rtilde_proposals(0.005, 0, 50, 1), 0.005)
})

test_that("negative binomial samplers adapt only at complete burn-in batches", {
  for (sampler in list(gibbs_SFM_neg_binomial, gibbs_SFM_shifted_neg_binomial)) {
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$draw_kap_mu = function(x, kappa, mu, ...) list(kappa = 0, mu = mu)
    sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
    for (burnin in c(0, 49, 100, 125)) {
      scales = numeric(0)
      batches = list()
      sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, mu_max) {
        scales <<- c(scales, std)
        iteration = length(scales)
        accepted = as.numeric(iteration <= 50 || iteration > 100)
        c(d_mu, d_rtilde, accepted, std)
      }
      sampler_env$adapt_rtilde_proposals = function(std, accepted, proposed, batch) {
        batches[[length(batches) + 1L]] <<- c(accepted, proposed, batch)
        adapt_rtilde_proposals(std, accepted, proposed, batch)
      }
      sampler(c(2, 3, 4), 1, 180, burnin, list(a0 = 1, A0 = 200, sig = 0.1), FALSE)
      expected = rep(0.1, 180)
      if (burnin >= 100) expected[51:100] = 0.1 * exp(0.1)
      expect_equal(scales, expected)
      expect_equal(batches, if (burnin >= 100) list(c(50, 50, 1), c(0, 50, 2)) else list())
    }
  }
})

test_that("rtilde batches exclude empty components and count only actual proposals", {
  for (sampler in list(gibbs_SFM_neg_binomial, gibbs_SFM_shifted_neg_binomial)) {
    sampler_env = new.env(parent = environment(sampler))
    environment(sampler) = sampler_env
    sampler_env$kmeans = function(y, ...) {
      list(cluster = rep(1L, length(y)), centers = matrix(c(2, 3, 4), ncol = 1))
    }
    iteration = 0L
    scales = numeric(0)
    batches = list()
    sampler_env$rdirichlet = function(n, alpha) {
      iteration <<- iteration + 1L
      matrix(if (iteration < 10) c(1, 0, 0) else c(0, 1, 0), nrow = n)
    }
    sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
    sampler_env$draw_kap_mu = function(x, kappa, mu, ...) list(kappa = 0, mu = mu)
    sampler_env$sample.int = function(...) 1L
    sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, mu_max) {
      scales <<- c(scales, std)
      c(d_mu, d_rtilde, as.numeric(length(Y) > 0 && d_mu == 2), std)
    }
    sampler_env$adapt_rtilde_proposals = function(std, accepted, proposed, batch) {
      batches[[length(batches) + 1L]] <<- list(accepted = accepted, proposed = proposed)
      adapt_rtilde_proposals(std, accepted, proposed, batch)
    }
    sampler(0:5, 3, 60, 50, list(a0 = 1, A0 = 200, sig = 0.1), FALSE)
    expect_equal(batches, list(list(accepted = c(10, 0, 0), proposed = c(10, 40, 0))))
    scales = matrix(scales, ncol = 3, byrow = TRUE)
    expect_equal(scales[1:50, ], matrix(0.1, nrow = 50, ncol = 3))
    expect_equal(scales[51:60, ], matrix(0.1 * exp(c(0.1, -0.1, 0)),
                                        nrow = 10, ncol = 3, byrow = TRUE))
  }
})

test_that("shifted negative binomial sorting preserves proposal scales and batch histories", {
  sampler = gibbs_SFM_shifted_neg_binomial
  sampler_env = new.env(parent = environment(sampler))
  environment(sampler) = sampler_env
  sampler_env$kmeans = function(y, ...) {
    list(cluster = c(1L, 1L, 2L, 2L), centers = matrix(c(2, 4), ncol = 1))
  }
  sampler_env$rdirichlet = function(n, alpha) matrix(c(0.5, 0.5), nrow = n)
  sampler_env$draw_e0 = function(e0, ...) list(e0, 0)
  allocation = 0L
  sampler_env$rmultinom = function(n, size, prob) {
    allocation <<- allocation %% 4L + 1L
    if (allocation <= 2L) c(1, 0) else c(0, 1)
  }
  shift_calls = 0L
  records = list()
  sampler_env$draw_kap_mu = function(x, kappa, mu, ...) {
    shift_calls <<- shift_calls + 1L
    iteration = ceiling(shift_calls / 2)
    reversed = iteration >= 25 && iteration < 75
    shift = if (mu == 2) as.numeric(reversed) else as.numeric(!reversed)
    list(kappa = shift, mu = mu)
  }
  sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, mu_max) {
    records[[length(records) + 1L]] <<- c(iteration = ceiling(shift_calls / 2),
                                        identity = d_mu, scale = std)
    c(d_mu, d_rtilde, as.numeric(d_mu == 2), std)
  }
  sampler(c(10, 11, 20, 21), 2, 160, 100, list(a0 = 1, A0 = 200, sig = 0.1), FALSE)
  records = do.call(rbind, records)
  direction = ifelse(records[, "identity"] == 2, 1, -1)
  adjustments = as.numeric(records[, "iteration"] > 50) + as.numeric(records[, "iteration"] > 100)
  expect_equal(unname(records[, "scale"]), 0.1 * exp(direction * adjustments * 0.1))
})

test_that("negative binomial mean draws use the fixed full-data bound", {
  update = rneg_binom
  update_env = new.env(parent = environment(update))
  environment(update) = update_env
  update_env$rgamma = function(n, ...) rep(1, n)
  update_env$rnorm = function(...) -1
  update_env$runif = function(...) 0.5

  zeros = update(c(0, 0), 0.1, 2, 0.2, mu_max = 10)
  small_counts = update(c(0, 1), 0.1, 2, 0.2, mu_max = 10)
  log_probability = pgamma(1 / 10, shape = 7, scale = 1 / 8,
                           lower.tail = FALSE, log.p = TRUE) + log(0.5)
  expected_mean = 1 / qgamma(log_probability, shape = 7, scale = 1 / 8,
                             lower.tail = FALSE, log.p = TRUE)
  expect_equal(zeros[1], expected_mean)
  expect_equal(small_counts[1], expected_mean)
  expect_gt(expected_mean, 1)
  expect_lt(expected_mean, 10)

  update_env$rgamma = function(n, ...) rep(1e5, n)
  result = update(c(0, 0), 0.1, 2, 0.2, mu_max = 10)
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

  result = update(numeric(0), 0.1, 2, 0.2, mu_max = 10)
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
    sampler_env$rneg_binom = function(Y, std, d_mu, d_rtilde, mu_max) {
      component_counts <<- c(component_counts, length(Y))
      prior_bounds <<- c(prior_bounds, mu_max)
      initial_means <<- c(initial_means, d_mu)
      rneg_binom(Y, std, d_mu, d_rtilde, mu_max = mu_max)
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