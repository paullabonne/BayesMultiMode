test_that("discrete_MF function returns expected results with dist = shifted_poisson and flat modes", {
  set.seed(1)
  lambda = c(0.1,1)
  kappa = c(10,0)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda, kappa = kappa)
  pars_names = c("eta", "lambda", "kappa")
  dist = "shifted_poisson"
  
  data = c(rpois(p[1]*1e3, lambda[1]) + kappa[1],
           rpois(p[2]*1e3, lambda[2]) + kappa[2])
  
  modes = discrete_MF(params, data = data, pars_names = pars_names, dist = dist)
  expect_equal(modes[1] == 0,  TRUE)
  expect_equal(modes[2] == 1,  TRUE)
  expect_equal(modes[3] == 10,  TRUE)
})

test_that("discrete_MF function returns expected results with dist = poisson", {
  set.seed(1)
  lambda = c(0.1,10)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda)
  pars_names = c("eta", "lambda")
  dist = "poisson"
  
  data = c(rpois(p[1]*1e3, lambda[1]),
           rpois(p[2]*1e3, lambda[2]))
  
  modes = discrete_MF(params, data = data, pars_names = pars_names, dist = dist)
  expect_equal(modes[1] == 0,  TRUE)
  expect_equal(modes[2] == 9,  TRUE)
})

test_that("discrete_MF function returns expected results with arbitrary function", {
  set.seed(1)
  mu = c(20,5)
  size = c(20,0.5)
  p = c(0.5,0.5)
  params = c(eta = p, mu = mu, size = size)
  pars_names = c("eta", "mu", "size")

  data = c(rnbinom(p[1]*1e3, mu = mu[1], size = size[1]),
           rnbinom(p[2]*1e3, mu = mu[2], size = size[2]))

  pmf_func <- function(x, pars) {
    dnbinom(x, mu = pars["mu"], size = pars["size"])
  }
  
  modes = discrete_MF(params, data = data, pars_names = pars_names, pmf_func = pmf_func)
  
  expect_equal(modes[1] == 0,  TRUE)
  expect_equal(modes[2] == 18,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
})