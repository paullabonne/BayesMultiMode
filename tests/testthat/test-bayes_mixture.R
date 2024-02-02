test_that("bayes_mixture returns expected results for the normal distribution", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.8,0.2)
  dist_type = "continuous"
  
  data = c(rnorm(p[1]*100, mu[1], sigma[1]), rnorm(p[2]*100, mu[2], sigma[2]))
  fit = c(eta = p, mu = mu, sigma = sigma)
  fit = rbind(fit, fit)

  BM = bayes_mixture(fit, data, burnin = 1, dist = "normal", pdf_func = NULL, dist_type = dist_type)
  expect_s3_class(BM, "bayes_mixture")
})

test_that("bayes_mixture returns expected results for the student distribution", {
  set.seed(1)
  mu = c(0.5,6)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.8,0.2)
  params = c(eta = p, mu = mu, sigma = sigma, nu = nu)
  dist_type = "continuous"
  
  data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
           sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
  
  fit = c(eta = p, mu = mu, sigma = sigma, nu = nu, xi = c(0,0))
  fit = rbind(fit, fit)
  
  pdf_func = function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }

  BM = bayes_mixture(
    fit,
    data,
    burnin = 1,
    pdf_func = pdf_func,
    dist_type = dist_type,
    loc = "mu"
  )
  expect_s3_class(BM, "bayes_mixture")
})