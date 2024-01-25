test_that("mix_mode() returns expected results with dist = shifted_poisson and flat modes", {
  set.seed(1)
  lambda = c(0.1,1)
  kappa = c(10,0)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda, kappa = kappa)
  dist = "shifted_poisson"
  
  data = c(rpois(p[1]*1e3, lambda[1]) + kappa[1],
           rpois(p[2]*1e3, lambda[2]) + kappa[2])
  
  mix = new_Mixture(params, data = data, dist = dist)
  modes = mix_mode(mix)
  
  # summary
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 1,  TRUE)
  expect_equal(modes$mode_estimates[3] == 10,  TRUE)
})

test_that("discrete_MF function returns expected results with dist = poisson", {
  set.seed(1)
  lambda = c(0.1,10)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda)
  dist = "poisson"
  
  data = c(rpois(p[1]*1e3, lambda[1]),
           rpois(p[2]*1e3, lambda[2]))
  
  mix = new_Mixture(params, data = data, dist = dist)
  modes = discrete_MF(mix)
  
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 9,  TRUE)
})

test_that("mix_mode() function returns expected results with arbitrary function", {
  set.seed(1)
  mu = c(20,5)
  size = c(20,0.5)
  p = c(0.5,0.5)
  params = c(eta = p, mu = mu, size = size)

  data = c(rnbinom(p[1]*1e3, mu = mu[1], size = size[1]),
           rnbinom(p[2]*1e3, mu = mu[2], size = size[2]))

  pmf_func <- function(x, pars) {
    dnbinom(x, mu = pars["mu"], size = pars["size"])
  }
  
  mix = new_Mixture(params, data = data, pdf_func = pmf_func,
                    dist_type = "discrete")
  modes = mix_mode(mix)
  
  # plot
  # expect_snapshot(plot(mix, from = 0, to = 50))
  # expect_snapshot(plot(modes, from = 0, to = 50))
  
  # summary
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 18,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
})

test_that("mix_mode() function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
  dist = "skew_normal"
  
  mix = new_Mixture(params, dist = dist)
  modes = mix_mode(mix)
  
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
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
  
  mix = new_Mixture(params, pdf_func = pdf_func, dist_type = "continuous")
  modes = mix_mode(mix)
  
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
})

test_that("mix_mode() function returns expected results", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.2,0.2)
  
  params = c(eta = p, mu = mu, sigma = sigma)
  
  mix = new_Mixture(params, dist = "normal")
  modes = mix_mode(mix)
  
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(round(modes$mode_estimates),  mu)
})