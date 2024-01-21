test_that("MEM function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
  dist = "skew_normal"

  mix = new_Mixture(params, dist = dist)
  modes = MEM(mix)
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
})

test_that("MEM function returns expected results with an arbitrary function", {
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
  modes = MEM(mix)
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
})