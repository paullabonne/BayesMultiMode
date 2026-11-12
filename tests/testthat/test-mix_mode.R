test_that("discrete modes exclude unrelated points with the same height", {
  pars = matrix(1, nrow = 1, dimnames = list(NULL, "eta"))
  pmf_func = function(x, pars) c(4, 3, 1, 2, 3, 2)[x + 1] / 15

  for (type in c("all", "unique")) {
    expect_warning(modes <- discrete_MF(pars, pmf_func, c(0, 5), type), NA)
    expect_equal(modes, c(0, 4))
  }
})

test_that("a flat peak has one mode and all of its locations", {
  probabilities = c(0.02, 0.04, 0.07, 0.10, 0.15, 0.20, 0.20, 0.10, 0.06, 0.04, 0.02)
  pmf_func = function(x, pars) probabilities[x + 1]
  params = c(eta = c(0.5, 0.5), mu = c(1, 2))
  mix = mixture(params, range = c(0, 10), pdf_func = pmf_func, dist_type = "discrete")

  expect_warning(modes <- mix_mode(mix), "Some modes are flat")
  expect_equal(modes$mode_estimates, c(5, 6))
  expect_output(summary(modes), "Number of modes found: 1")
  expect_warning(unique_modes <- mix_mode(mix, type = "unique"), "Some modes are flat")
  expect_length(unique_modes$mode_estimates, 1)
  expect_output(summary(unique_modes), "Number of modes found: 1")

  fit = rbind(params, params, params, params)
  bayesmix = bayes_mixture(fit, data = 0:10, burnin = 1,
                          pdf_func = pmf_func, dist_type = "discrete")
  bayesmodes = suppressWarnings(bayes_mode(bayesmix))
  expect_equal(bayesmodes$p1, 1)
  expect_equal(unname(bayesmodes$p_mode_loc), rbind(c(5, 6), c(1, 1)))
})

test_that("experiment helpers keep plateau locations separate from peak counts", {
  for (script_name in c("run-validation.R", "run-comparison.R")) {
    script_path = test_path("..", "..", "JOE revision", "experiment", script_name)
    skip_if_not(file.exists(script_path), "Experiment scripts are not included in the package")
    expressions = parse(script_path)
    helper_assignment = vapply(expressions, function(expression) {
      is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        identical(expression[[2]], as.name("pmf_modes"))
    }, logical(1))
    expect_equal(sum(helper_assignment), 1)
    helper_environment = new.env(parent = baseenv())
    eval(expressions[[which(helper_assignment)]], helper_environment)
    find_modes = helper_environment$pmf_modes

    probabilities = c(0.02, 0.04, 0.07, 0.10, 0.15, 0.20, 0.20, 0.10, 0.06, 0.04, 0.02)
    peaks = find_modes(0:10, probabilities)
    expect_equal(peaks, list(c(5L, 6L)))
    expect_length(peaks, 1)
    draws = list(peaks, peaks)
    expect_equal(lengths(draws), c(1L, 1L))
    expect_equal(as.numeric(table(unlist(draws))) / length(draws), c(1, 1))
    expect_equal(find_modes(0:7, c(0.1, 0.4, 0.4, 0.1, 0.2, 0.5, 0.5, 0.1)),
                 list(c(1L, 2L), c(5L, 6L)))
    expect_equal(find_modes(0:2, c(0.4, 0.4, 0.2)), list(c(0L, 1L)))
    expect_equal(find_modes(0:3, c(0.1, 0.4, 0.1, 0.3)), list(1L, 3L))
    expect_identical(find_modes(0:2, c(0, 0, 0)), list())
  }
})

test_that("mix_mode() returns expected results with dist = shifted_poisson and flat modes", {
  set.seed(1)
  lambda = c(0.1,1)
  kappa = c(10,0)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda, kappa = kappa)
  dist = "shifted_poisson"
  
  mix = mixture(params, range = c(0, 50), dist = dist)
  modes = mix_mode(mix)
  
  # summary
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 1,  TRUE)
  expect_equal(modes$mode_estimates[3] == 10,  TRUE)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with dist = poisson", {
  set.seed(1)
  lambda = c(0.1,10)
  p = c(0.5,0.5)
  params = c(eta = p, lambda = lambda)
  dist = "poisson"
  
  mix = mixture(params, range = c(0, 50), dist = dist)
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 9,  TRUE)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with arbitrary function", {
  set.seed(1)
  mu = c(20,5)
  size = c(20,0.5)
  p = c(0.5,0.5)
  params = c(eta = p, mu = mu, size = size)


  pmf_func <- function(x, pars) {
    dnbinom(x, mu = pars["mu"], size = pars["size"])
  }
  
  mix = mixture(params, range = c(0, 50), pdf_func = pmf_func,
                    dist_type = "discrete")
  modes = mix_mode(mix)
  
  # plot
  # expect_snapshot(plot(mix, from = 0, to = 50))
  # expect_snapshot(plot(modes, from = 0, to = 50))
  
  # summary
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(modes$mode_estimates[1] == 0,  TRUE)
  expect_equal(modes$mode_estimates[2] == 18,  TRUE)
  # the two densities are far appart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with dist = skew_normal", {
  set.seed(1)
  xi = c(0,6)
  omega = c(1,2)
  alpha = c(0,0)
  p = c(0.8,0.2)
  params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
  dist = "skew_normal"
  
  mix = mixture(params, dist = dist, range = c(-5,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results with an arbitrary function", {
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
  
  mix = mixture(params, pdf_func = pdf_func,
                    dist_type = "continuous",
                    loc = "mu",
                range = c(-5,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(abs(sum(modes$mode_estimates-xi))<1,  TRUE)
  # the two densities are far apart so the modes should coincide with the location parameters
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})

test_that("mix_mode() function returns expected results", {
  set.seed(1)
  mu = c(0,5)
  sigma = c(1,2)
  p = c(0.8,0.2)
  
  params = c(eta = p, mu = mu, sigma = sigma)
  
  mix = mixture(params, dist = "normal", range = c(-4,10))
  modes = mix_mode(mix)
  
  expect_snapshot(modes$mode_estimates)
  
  expect_equal(round(modes$mode_estimates),  mu)
  
  skip_on_ci()
  expect_snapshot(summary(mix))
  expect_snapshot(summary(modes))
})