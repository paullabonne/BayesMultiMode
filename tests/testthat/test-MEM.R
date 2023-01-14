test_that("MEM function returns expected results with dist = student", {
  set.seed(1)
  mu = c(0.5,6)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.8,0.2)
  params = c(theta = p, mu = mu, sigma = sigma, nu = nu)
  pars_names = c("theta", "mu", "sigma", "nu")
  
  data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
           sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
  
  modes = MEM(params, "student", pars_names, data)
  expect_equal(abs(sum(modes-mu))<1,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
})

test_that("MEM function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(theta = p, mu = xi, sigma = omega, xi = alpha)
  pars_names = c("theta", "mu", "sigma", "xi")
  
  data = c(sn::rsn(p[1]*100, xi[1], omega[1], alpha = alpha[1]),
           sn::rsn(p[2]*100, xi[2], omega[2], alpha = alpha[2]))
  
  modes = MEM(params, "skew_normal", pars_names, data)
  expect_equal(abs(sum(modes-xi))<1,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
})