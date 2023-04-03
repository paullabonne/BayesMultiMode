
test_that("fixed_point function returns expected results", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.2,0.2)
  
  data = c(rnorm(p[1]*100, mu[1], sigma[1]), rnorm(p[2]*100, mu[2], sigma[2]))
  params = c(eta = p, mu = mu, sigma = sigma)
  pars_names = c("eta", "mu", "sigma")
  
  modes = fixed_point(params, data, pars_names)
  expect_equal(round(modes),  mu)
})