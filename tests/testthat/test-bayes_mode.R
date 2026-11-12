test_that("zero-inflated fits support full and positive-count modes and plots", {
  observations = c(rep(0, 20), rep(4:12, 5))
  for (distribution in c("poisson", "neg_binomial")) {
    set.seed(615)
    fit = bayes_fit(observations, 1, distribution, nb_iter = 60, burnin = 30,
                    print = FALSE, zero_inflated = TRUE)
    modes = bayes_mode(fit)
    expect_s3_class(modes, "bayes_mode")
    expect_true(all(modes$modes[, 1] == 0))
    expect_true(all(is.finite(as.matrix(modes$mix_density))))
    positive_modes = bayes_mode(fit, range = c(1, max(observations)))
    expect_true(all(positive_modes$modes[!is.na(positive_modes$modes)] >= 1))
    expect_equal(sum(positive_modes$p_nb_modes[2, ]), 1)
    for (graph in list(plot(fit, draws = 3), plot(modes), plot(modes, graphs = "loc"))) {
      expect_s3_class(graph, "ggplot")
      expect_no_error(ggplot2::ggplot_build(graph))
    }
    expect_no_error(plot(fit, draw = 1))
    expect_no_error(plot(modes, draw = 1))
    expect_output(summary(fit), "Summary of MCMC output")
    expect_output(summary(modes), "Posterior probability")
  }
})

test_that("the zero atom is an ordinary discrete mode in imported draws", {
  for (distribution in c("poisson", "neg_binomial")) {
    parameters = if (distribution == "poisson") {
      c(eta1 = 0.35, eta2 = 0.65, lambda1 = 0, lambda2 = 7.5)
    } else {
      c(eta1 = 0.35, eta2 = 0.65, mu1 = 0, mu2 = 8,
         rtilde1 = 0.25, rtilde2 = 0.1)
    }
    fit = bayes_mixture(rbind(parameters, parameters), 0:25, dist = distribution)
    modes = bayes_mode(fit)
    expect_equal(unname(modes$modes[, 1:2]), matrix(c(0, 0, 7, 7), nrow = 2))
    expect_true(all(is.na(modes$modes[, -(1:2), drop = FALSE])))
    expect_equal(unname(modes$p_nb_modes), matrix(c(2, 1), ncol = 1))
  }
})

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