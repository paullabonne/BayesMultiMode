test_that("bayes_mode works with external MCMC output", {
  set.seed(1)
  mu = c(0.5,6)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.8,0.2)
  params = c(eta = p, mu = mu, sigma = sigma, nu = nu)
  pars_names = c("eta", "mu", "sigma", "nu")
  dist_type = "continuous"
  
  data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
           sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
  
  fit = c(eta = p, mu = mu, sigma = sigma, nu = nu)
  fit = rbind(fit, fit)
  
  pdf_func = function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }
  
  bayesmix = new_BayesMixture(fit, data, K = 2, burnin = 1, pars_names = pars_names,
                                pdf_func = pdf_func, dist_type = dist_type)
  
  bayesmode = bayes_mode(bayesmix)
  
  expect_equal(abs(sum(bayesmode$modes-mu))<0.1,  TRUE)
})