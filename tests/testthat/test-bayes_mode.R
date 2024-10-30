test_that("bayes_mode works with external MCMC output", {
  set.seed(1)
  mu = c(0.5,6)
  sigma = c(1,2)
  nu = c(5,5)
  p = c(0.8,0.2)
  dist_type = "continuous"
  
  data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
           sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
  
  fit = c(eta = p, mu = mu, sigma = sigma, xi = c(0,0), nu = nu)
  fit = rbind(fit, fit)
  
  pdf_func = function(x, pars) {
    sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
  }
  
  bayesmix = bayes_mixture(fit,
                           data,
                           burnin = 1,
                           pdf_func = pdf_func,
                           dist_type = dist_type,
                           loc = "mu")
  
  bayesmode = bayes_mode(bayesmix)
  m = apply(bayesmode$modes,2,mean)
  m = m[order(m)]
  m = m[!is.na(m)]
  
  expect_equal(sum(abs(m-mu)<0.1),  2)
  
  skip_on_ci()
  expect_snapshot(summary(bayesmix))
  expect_snapshot(summary(bayesmode))
})

test_that("bayes_mode works with normal mixture", {
  set.seed(123)
  
  mu = c(-5,5)
  # retrieve galaxy data
  y = rnorm(200, mu)
  
  # estimation
  bayesmix = bayes_fit(data = y,
                       K = 2, #not many to run the example rapidly
                       dist = "normal",
                       nb_iter = 500, #not many to run the example rapidly
                       burnin = 100)
  
  # mode estimation
  bayesmode = bayes_mode(bayesmix)
  m = apply(bayesmode$modes,2,mean)
  m = m[order(m)]
  m = m[!is.na(m)]
  
  expect_equal(sum(abs(m-mu)<0.1),  2)
  
  skip_on_ci()
  expect_snapshot(summary(bayesmix))
  expect_snapshot(summary(bayesmode))
})

test_that("bayes_mode works with skew_normal mixture", {
  set.seed(123)
  
  mu = c(-5,5)
  # retrieve galaxy data
  y = rnorm(200, mu)
  
  # estimation
  bayesmix = bayes_fit(data = y,
                       K = 2, #not many to run the example rapidly
                       dist = "skew_normal",
                       nb_iter = 500, #not many to run the example rapidly
                       burnin = 100)
  
  # mode estimation
  bayesmode = bayes_mode(bayesmix)
  m = apply(bayesmode$modes,2,mean, na.omit = T)
  m = m[order(m)]
  m = m[!is.na(m)]
  
  skip_on_cran()
  expect_equal(sum(abs(m-mu)<0.5),  2)
  
  skip_on_ci()
  expect_snapshot(summary(bayesmix))
  expect_snapshot(summary(bayesmode))
})

test_that("bayes_mode works with shifted poisson mixture", {
  set.seed(123)
  
  mu = c(0,5)
  # data
  y = c(rpois(100, 1) + mu[1],
        rpois(100, 1) + mu[2])
  
  # estimation
  bayesmix = bayes_fit(data = y,
                       K = 2, #not many to run the example rapidly
                       dist = "shifted_poisson",
                       nb_iter = 500, #not many to run the example rapidly
                       burnin = 100)
  
  # mode estimation
  bayesmode = bayes_mode(bayesmix)
  
  # check the plots and summary
  # plot(bayesmix)
  
  # check the modes
  m = apply(bayesmode$modes,2,median, na.rm=T)
  m = m[order(m)]
  m = m[!is.na(m)]
  
  expect_equal(sum(abs(m-mu)<0.5),  2)
  
  skip_on_ci()
  expect_snapshot(summary(bayesmix))
  expect_snapshot(summary(bayesmode))
})

test_that("bayes_mode works with poisson mixture", {
  set.seed(123)
  
  # retrieve galaxy data
  y = c(rpois(100, 0.5),
        rpois(100, 10))
  
  # estimation
  bayesmix = bayes_fit(data = y,
                       K = 2, #not many to run the example rapidly
                       dist = "poisson",
                       nb_iter = 500, #not many to run the example rapidly
                       burnin = 100)
  
  # mode estimation
  bayesmode = bayes_mode(bayesmix)
  
  # check the plots and summary
  # plot(bayesmix)
  
  # check the modes
  m = apply(bayesmode$modes,2,median, na.rm=T)
  m = m[order(m)]
  m = m[!is.na(m)]
  
  expect_equal(sum(abs(m-c(0,10))<2),  2)
  
  skip_on_ci()
  expect_snapshot(summary(bayesmix))
  expect_snapshot(summary(bayesmode))
})