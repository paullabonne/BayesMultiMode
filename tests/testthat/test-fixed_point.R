
test_that("fixed_point function returns expected results", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.2,0.2)
  
  params = c(eta = p, mu = mu, sigma = sigma)
  pars_names = c("eta", "mu", "sigma")
  
  mix = new_Mixture(params, pars_names = pars_names)
  modes = fixed_point(mix)
  expect_equal(round(modes$mode_estimates),  mu)
})