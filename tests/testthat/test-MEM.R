test_that("MEM function returns expected results", {
  set.seed(1)
  mu = c(0,10)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.2,0.2)
  
  data = c(sn::rst(p[1]*100, mu[1], sigma[1], nu = nu[1]), sn::rst(p[2]*100, mu[2], sigma[2], nu = nu[2]))
  params = c(theta = p, mu = mu, sigma = sigma, nu = nu)
  
  modes = MEM(params, "student", data)
})