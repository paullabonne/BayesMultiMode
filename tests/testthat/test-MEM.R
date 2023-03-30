test_that("MEM function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
  pars_names = c("eta", "xi", "omega", "alpha")
  dist = "skew_normal"
  
  data = c(sn::rsn(p[1]*1e3, xi[1], omega[1], alpha[1]),
           sn::rsn(p[2]*1e3, xi[2], omega[2], alpha[2]))
  
  modes = MEM(params, data = data, pars_names = pars_names, dist = dist)
  expect_equal(abs(sum(modes-xi))<1,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
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
  pars_names = c("eta", "mu", "sigma", "xi", "nu")
  
  pdf_func <- function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }
  
  data = c(sn::rst(p[1]*1e3, xi[1], omega[1], alpha = alpha[1], nu = nu[1]),
           sn::rst(p[2]*1e3, xi[2], omega[2], alpha = alpha[2], nu = nu[2]))
  
  modes = MEM(params, pars_names = pars_names, data = data, pdf_func = pdf_func)
  expect_equal(abs(sum(modes-xi))<1,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
})