test_that("mix_mode() returns expected results with dist = shifted_poisson and flat modes", {
  set.seed(1)
  lambda = c(0.1,1)
  kappa = c(10,0)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda, kappa = kappa)
  dist = "shifted_poisson"
  
  mix = mixture(params, range = c(0, 50), dist = dist)
  modes = mix_mode(mix)
  
  # summary
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 1,  TRUE)
  expect_equal(modes$mode_estimates[3] == 10,  TRUE)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with dist = poisson", {
  set.seed(1)
  lambda = c(0.1,10)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda)
  dist = "poisson"
  
  mix = mixture(params, range = c(0, 50), dist = dist)
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 9,  TRUE)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with arbitrary function", {
  set.seed(1)
  mu = c(20,5)
  size = c(20,0.5)
  p = c(0.5,0.5)
  params = c(eta = p, mu = mu, size = size)


  pmf_func <- function(x, pars) {
    dnbinom(x, mu = pars["mu"], size = pars["size"])
  }
  
  mix = mixture(params, range = c(0, 50), pdf_func = pmf_func,
                    dist_type = "discrete")
  modes = mix_mode(mix)
  
  # plot
  # expect_snapshot(plot(mix, from = 0, to = 50))
  # expect_snapshot(plot(modes, from = 0, to = 50))
  
  # summary
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 18,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
  dist = "skew_normal"
  
  mix = mixture(params, dist = dist, range = c(-5,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with an arbitrary function", {
  # example with the skew-t of the sn package
  
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  nu = c(3,100)
  p = c(0.8,0.2)
  params = c(eta = p, mu = xi, sigma = omega, xi = alpha, nu = nu)
  
  pdf_func <- function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }
  
  mix = mixture(params, pdf_func = pdf_func,
                    dist_type = "continuous",
                    loc = "mu",
                range = c(-5,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.8,0.2)
  
  params = c(eta = p, mu = mu, sigma = sigma)
  
  mix = mixture(params, dist = "normal", range = c(-4,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(round(modes$mode_estimates),  mu)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})