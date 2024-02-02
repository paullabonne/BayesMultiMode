test_that("mixture() returns expected error when dist and pdf_func are not given", {
  
  params = c(eta = c(0.2,0.2),
             mu = c(0,5),
             sigma = c(1,2))
  
  expect_error(mixture(params),
               "one of dist or pdf_func must be specified")
})

test_that("mixture() returns expected error when dist and parameters do not match (skew normal)", {
  
  params = c(eta = c(0.8,0.2),
             mu = c(0,5),
             omega = c(1,2),
             alpha = c(0.5,0.1))
  
  expect_error(mixture(params, dist = "skew_normal", dist_type = "continuous"),
               "variable names in pars should be eta, xi, omega and alpha when dist = skew_normal")
})

test_that("mixture() returns expected error when pdf_func and mcmc parameters do not match", {
  set.seed(1)
  mu = c(0.5,6)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.8,0.2)
  
  # name of scale set to omega instead of sigma
  fit = c(eta = p, mu = mu, omega = sigma, nu = nu, xi = c(0,0))

  pdf_func <- function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }
  
  expect_error(mixture(fit, pdf_func = pdf_func,
                           dist_type = "continuous",
                           loc = "mu"),
               "running pdf_func for the first component returns NA")

})