test_that("telescoping conditionals agree with a labelled partition enumeration", {
  probability <- function(counts, total, shape) {
    labels <- length(counts)
    assignments <- prod(seq(total, length.out = labels, by = -1))
    numerator <- prod(vapply(counts, function(size) prod(shape + seq_len(size) - 1), 0.0))
    assignments * numerator / prod(total * shape + seq_len(sum(counts)) - 1)
  }
  totals <- 2:5
  concentration <- 1.7
  normalize <- function(mass) mass / sum(mass)
  for (model in c("static", "dynamic")) {
    shape <- function(total) if (model == "static") concentration else concentration / total
    for (prior in c("uniform", "bnb")) {
      reference <- vapply(totals, function(total) {
        prior_mass <- if (prior == "uniform") 1 else integrate(function(prob)
          stats::dnbinom(total - 1, 1, prob) * stats::dbeta(prob, 4, 3), 0, 1)$value
        prior_mass * probability(c(3, 1), total, shape(total))
      }, 0.0)
      actual <- function(counts) normalize(exp(tele_log_count(totals, counts, concentration, model, prior)))
      expect_equal(actual(c(3, 1)), normalize(reference))
      expect_identical(isTRUE(all.equal(actual(c(3, 1)), actual(c(2, 2)))), model == "static")
    }
    values <- c(0.7, concentration)
    reference <- vapply(values, function(value) {
      shape <- if (model == "static") value else value / 4
      density <- if (model == "static") stats::dgamma(value, shape = 1, rate = 20) else stats::df(value, 6, 3)
      probability(c(3, 1), 4, shape) * density * value
    }, 0.0)
    target <- vapply(log(values), tele_log_concentration, 0.0, total = 4, counts = c(3, 1), model = model)
    expect_equal(diff(target), diff(log(reference)))
  }
})

test_that("telescoping follows the concentration defaults of the telescope prior helpers", {
  draws <- function(model, ...) {
    set.seed(206)
    fit <- bayes_fit(c(1, 2, 3), 1, "normal", nb_iter = 50, burnin = 0, print = FALSE,
                     sampler = "telescoping", weight_model = model, ...)
    unname(fit$state_all[, "concentration"])
  }
  expect_identical(draws("dynamic"), rep(1, 50))
  dynamic <- draws("dynamic", concentration = "F_6_3")
  expect_gt(length(unique(dynamic)), 1L)
  expect_identical(dynamic, draws("dynamic", concentration = NULL))
  static <- draws("static")
  expect_gt(length(unique(static)), 1L)
  expect_identical(static, draws("static", concentration = "G_1_20"))
  expect_error(draws("static", concentration = "F_6_3"), "concentration")
})

test_that("occupied labels carry their allocations and component state", {
  allocations <- diag(4)[c(4, 2, 4), ]
  pars <- cbind(mu = c(10, 20, 30, 40), scale = 1:4)
  occupied <- tele_occupied(allocations, pars)
  expect_equal(occupied$counts, c(1, 2))
  expect_equal(occupied$pars, pars[c(2, 4), ])
  expect_equal(occupied$allocations, allocations[, c(2, 4)])
  set.seed(71)
  values <- 0
  family <- tele_normal_family(values, check_priors(list(B0 = 1), "normal", values))
  fit <- gibbs_telescoping(values, 5, 5, 200, family, "static", "uniform", 1, FALSE)
  expect_equal(range(fit$state[, "K"]), c(1, 5))
  expect_true(any(fit$state[, "K"] > fit$state[, "Kplus"]))
  expect_equal(rowSums(!is.na(fit$draws)), 3 * fit$state[, "K"])
  expect_equal(rowSums(fit$draws[, 1:5], na.rm = TRUE), rep(1, 200))
})

test_that("public draws preserve the density, likelihood and burn-in alignment", {
  set.seed(64)
  values <- c(-0.2, 0, 0.2)
  cap <- 5
  fit <- bayes_fit(values, cap, "normal", nb_iter = 70, burnin = 10, print = FALSE,
                   sampler = "telescoping", K_init = cap, concentration = "F_6_3")
  expect_equal(fit$K, cap)
  expect_identical(fit$state, fit$state_all[11:70, ])
  expect_identical(fit$mcmc, fit$mcmc_all[11:70, ])
  expect_equal(fit$cap_mass, mean(fit$state[, "K"] == cap))
  expect_true(any(fit$state[, "K"] > fit$state[, "Kplus"]))
  expect_gt(length(unique(fit$state[, "concentration"])), 1)
  reference <- vapply(seq_len(nrow(fit$mcmc_all)), function(iteration) {
    total <- fit$state_all[iteration, "K"]
    pars <- matrix(fit$mcmc_all[iteration, ], nrow = cap)
    active <- seq_len(total)
    stopifnot(all(is.na(pars[-active, , drop = FALSE])), all(is.finite(pars[active, ])))
    Reduce(`+`, lapply(active, function(component) {
      pars[component, 1] * dnorm(values, pars[component, 2], pars[component, 3])
    }))
  }, numeric(length(values)))
  expect_equal(colSums(log(reference)), fit$loglik)
  reconstructed <- apply(fit$mcmc, 1, dmix, x = values, pars_names = fit$pars_names, pdf_func = fit$pdf_func)
  expect_equal(reconstructed, reference[, 11:70])
  expect_output(summary(fit), "Posterior mass at the cap")
  trace <- suppressMessages(bayes_trace(fit, c("K", "Kplus", "concentration"), with_burnin = TRUE))
  expect_s3_class(trace, "ggplot")
  expect_equal(nrow(trace$data), 3 * nrow(fit$state_all))
})

test_that("all telescope families preserve their public density with either weight model", {
  fields <- list(normal = c("mu", "sigma"), skew_normal = c("xi", "omega", "alpha"),
                 poisson = "lambda", shifted_poisson = c("kappa", "lambda"),
                 neg_binomial = c("mu", "rtilde"), shifted_neg_binomial = c("kappa", "mu", "rtilde"))
  cap <- 5L
  for (dist in names(fields)) {
    values <- if (dist %in% c("normal", "skew_normal")) c(-2, -0.5, 0, 0.5, 3) else c(0, 0, 0, 1, 4)
    expect_error(bayes_fit(values, cap, dist, sampler = "telescoping", zero_inflated = TRUE), "zero inflation.*sfm")
    for (model in c("dynamic", "static")) {
      set.seed(389)
      fit <- bayes_fit(values, cap, dist, nb_iter = 60, burnin = 10, print = FALSE,
                       sampler = "telescoping", weight_model = model, K_prior = "uniform", K_init = cap)
      total <- fit$state_all[, "K"]
      occupied <- fit$state_all[, "Kplus"]
      expect_true(all(occupied >= 1 & occupied <= total & total <= cap))
      expect_true(any(occupied < total))
      expect_identical(fit$pars_names, c("eta", fields[[dist]]))
      slots <- rep(seq_len(cap), times = length(fit$pars_names))
      expect_identical(is.na(unname(fit$mcmc_all)), outer(total, slots, `<`))
      expect_equal(rowSums(fit$mcmc_all[, seq_len(cap)], na.rm = TRUE), rep(1, length(total)))
      reference <- apply(fit$mcmc_all, 1, function(draw) {
        sum(log(dmix(values, draw, pars_names = fit$pars_names, pdf_func = fit$pdf_func)))
      })
      expect_equal(fit$loglik, reference)
    }
  }
})

test_that("one-component Poisson rates follow the conjugate posterior", {
  values <- c(0, 0, 0, 1, 4)
  priors <- list(l0 = 2, L0 = 3)
  set.seed(241)
  fit <- bayes_fit(values, 1, "poisson", priors, nb_iter = 1000, burnin = 0,
                   print = FALSE, sampler = "telescoping")
  shape <- priors$l0 + sum(values)
  rate <- priors$L0 + length(values)
  draws <- fit$mcmc[, "lambda1"]
  expect_lt(abs(mean(draws) - shape / rate), 5 * sqrt(shape / length(draws)) / rate)
  expect_lt(abs(mean(draws <= qgamma(0.75, shape, rate)) - 0.75), 5 * sqrt(0.75 * 0.25 / length(draws)))
})

test_that("empty skew-normal components follow their priors and use residual variance for the shared rate", {
  priors <- check_priors(list(b0 = -3, D_xi = 4, D_psi = 0.25), "skew_normal", 0:4)
  family <- tele_skew_normal_family(0:4, priors)
  set.seed(423)
  pars <- family$empty(2000, 3)
  expect_lt(abs(mean(pars[, "xi"]) - priors$b0), 5 * sqrt(priors$D_xi / nrow(pars)))
  expect_lt(abs(mean(pars[, "psi"]^2) - priors$D_psi), 5 * priors$D_psi * sqrt(2 / nrow(pars)))
  expect_lt(abs(mean(1 / pars[, "variance"]) - priors$c0 / 3), 5 * sqrt(priors$c0 / nrow(pars)) / 3)
  expect_equal(pars[, "omega"]^2, pars[, "variance"] + pars[, "psi"]^2)
  expect_equal(pars[, "alpha"], pars[, "psi"] / sqrt(pars[, "variance"]))
  occupied <- pars[1:3, ]
  set.seed(66)
  expected <- rgamma(1, priors$g0 + 3 * priors$c0, priors$G0 + sum(1 / occupied[, "variance"]))
  set.seed(66)
  expect_equal(family$shared(occupied), expected)
})

test_that("telescoping default priors remain proper for degenerate data", {
  expect_equal(check_priors(list(), "poisson", c(0, 0, 4), telescoping = TRUE)$L0, 1.1 / 4)
  expect_equal(check_priors(list(), "poisson", 0, telescoping = TRUE)$L0, 1.1)
  for (dist in c("poisson", "shifted_poisson", "skew_normal")) {
    set.seed(390)
    fit <- bayes_fit(0, 4, dist, nb_iter = 15, burnin = 0, print = FALSE, sampler = "telescoping", K_init = 4)
    expect_true(all(is.finite(fit$loglik)))
  }
  expect_error(bayes_fit(0, 1, "neg_binomial", sampler = "telescoping"), "full-data maximum")
})

test_that("one-component draws pass through mode inference and plotting", {
  set.seed(10)
  fit <- bayes_fit(c(1, 2, 3), 1, "normal", nb_iter = 12, burnin = 2, print = FALSE,
                   sampler = "telescoping", weight_model = "static", concentration = 1)
  modes <- bayes_mode(fit, inside_range = FALSE)
  expect_equal(modes$p1, 1)
  expect_equal(as.numeric(modes$modes), unname(fit$mcmc[, "mu1"]))
  expect_true(all(fit$state == 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, draw = 1))
  expect_no_error(plot(modes, draw = 1))
  expect_output(summary(modes), "multimodality")
  padded <- c(eta1 = 0.8, eta2 = 0.2, eta3 = NA, mu1 = 0, mu2 = 12, mu3 = NA,
              sigma1 = 1, sigma2 = 1, sigma3 = NA)
  mix <- mixture(padded, dist = "normal", range = c(-5, 15))
  expect_equal(mix$K, 2)
  expect_equal(unname(mix_mode(mix, tol_mixp = 0.9)$mode_estimates), 0)
})

test_that("telescoping controls are checked without changing seeded default SFM", {
  values <- c(-2, -1, 0, 1, 2, 3)
  set.seed(12)
  fit <- bayes_fit(values, 2, "normal", list(), 20, 10, FALSE, FALSE)
  expect_equal(unname(c(fit$mcmc_all[20, ], fit$loglik[20])),
               c(0.658603976534, 0.341396023466, -1.11754146012, 4.10743629172,
                 1.16904290244, 2.37435915408, -13.8476242980))
  set.seed(12)
  expect_identical(bayes_fit(values, 2, "normal", list(), 20, 10, FALSE, FALSE, "sfm")$mcmc_all, fit$mcmc_all)
  run <- function(dist = "normal", ...) bayes_fit(1, 5, dist, nb_iter = 4, burnin = 0,
                                                  print = FALSE, sampler = "telescoping", ...)
  expect_no_error(run(K_init = 5))
  expect_error(run("unknown"), "Unsupported distribution")
  expect_error(run(zero_inflated = TRUE), "zero inflation")
  expect_error(run(K_init = 6), "K_init")
  expect_error(run(concentration = 0), "concentration")
})

test_that("shifted NB fits retain empty components and use the full-data mean bound", {
  set.seed(201)
  values <- c(0, 1, 7, 8)
  fit <- bayes_fit(values, 8, "shifted_neg_binomial", nb_iter = 300, burnin = 100,
                   print = FALSE, sampler = "telescoping", K_init = 8)
  expect_true(any(fit$state[, "K"] > fit$state[, "Kplus"]))
  expect_equal(fit$pars_names, c("eta", "kappa", "mu", "rtilde"))
  reconstructed <- vapply(seq_len(nrow(fit$state_all)), function(iteration) {
    total <- fit$state_all[iteration, "K"]
    pars <- matrix(fit$mcmc_all[iteration, ], nrow = 8)[seq_len(total), , drop = FALSE]
    stopifnot(all(pars[, 3] > 0 & pars[, 3] <= max(values)), all(pars[, 4] > 0 & pars[, 4] < 0.5))
    mass <- Reduce(`+`, lapply(seq_len(total), function(component) {
      pars[component, 1] * dnbinom(values - pars[component, 2],
        size = 1 / pars[component, 4] - 1, mu = pars[component, 3])
    }))
    sum(log(mass))
  }, 0.0)
  expect_equal(reconstructed, fit$loglik)
  expect_error(bayes_fit(0, 2, "shifted_neg_binomial", sampler = "telescoping"), "full-data maximum")
  expect_error(bayes_fit(values, 2, "shifted_neg_binomial", sampler = "telescoping", zero_inflated = TRUE), "zero inflation")
})

test_that("NB proposal histories follow occupied rows and freeze after burn-in", {
  set.seed(292)
  family <- tele_nb_family(c(0, 8), list(sig = 0.13), 50, shifted = TRUE)
  pars <- family$empty(3, NULL)
  expect_equal(unname(pars[, "scale"]), rep(0.13, 3))
  pars[, "scale"] <- c(0.2, 0.25, 0.3)
  pars[, "accepted"] <- c(40, 20, 1)
  pars[, "proposed"] <- 49
  occupied <- tele_occupied(diag(3)[c(3, 1), ], pars)$pars
  for (component in seq_len(nrow(occupied))) {
    current <- occupied[component, ]
    adapted <- family$occupied(8, current, NULL, 50)
    expected <- current[["scale"]] * exp(if (component == 1) 0.1 else -0.1)
    expect_equal(adapted[["scale"]], expected)
    expect_equal(unname(adapted[c("accepted", "proposed")]), c(0, 0))
    frozen <- family$occupied(8, adapted, NULL, 100)
    expect_identical(frozen[c("scale", "accepted", "proposed")], adapted[c("scale", "accepted", "proposed")])
  }
  expect_equal(unname(family$empty(1, NULL)[1, c("scale", "accepted", "proposed")]), c(0.13, 0, 0))
})