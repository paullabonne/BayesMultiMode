simulation_experiments <- list(
  normal_unimodal = list(dist = "normal", location = c(-0.5, 0.5), scale = 1,
                         grid = seq(-12, 12, by = 0.025), modes = 1L),
  normal_bimodal = list(dist = "normal", location = c(-3, 3), scale = 0.7,
                        grid = seq(-12, 12, by = 0.025), modes = 2L),
  skew_normal_unimodal = list(dist = "skew_normal", location = c(-0.25, 0.25), scale = 1,
                              shape = c(3, 3), grid = seq(-12, 12, by = 0.025), modes = 1L),
  skew_normal_bimodal = list(dist = "skew_normal", location = c(-4, 4), scale = 0.8,
                             shape = c(4, -4), grid = seq(-12, 12, by = 0.025), modes = 2L),
  poisson_unimodal = list(dist = "poisson", rate = c(4.5, 6.5), grid = 0:200, modes = 1L),
  poisson_bimodal = list(dist = "poisson", rate = c(0.5, 10.5), grid = 0:200, modes = 2L),
  shifted_poisson_unimodal = list(dist = "shifted_poisson", rate = c(20.5, 10.5),
                                  shift = c(0, 10), grid = 0:200, modes = 1L),
  shifted_poisson_bimodal = list(dist = "shifted_poisson", rate = c(2.5, 2.5),
                                 shift = c(0, 10), grid = 0:200, modes = 2L),
  neg_binomial_unimodal = list(dist = "neg_binomial", mean = c(2, 12), size = 4,
                               grid = 0:200, modes = 1L),
  neg_binomial_bimodal = list(dist = "neg_binomial", mean = c(1.5, 24), size = 6,
                              grid = 0:200, modes = 2L),
  shifted_neg_binomial_unimodal = list(dist = "shifted_neg_binomial", mean = c(16, 6),
                                       size = 6, shift = c(0, 10), grid = 0:200, modes = 1L),
  shifted_neg_binomial_bimodal = list(dist = "shifted_neg_binomial", mean = c(2, 7.5),
                                      size = 4, shift = c(0, 10), grid = 0:200, modes = 2L),
  shifted_poisson_comparison_c4 = list(dist = "shifted_poisson", rate = c(5, 10),
                                       shift = c(10, 30), grid = 0:200, modes = 2L,
                                       fit = list(count = 500L, K = 6L,
                                                  nb_iter = 5000L, burnin = 2500L)),
  shifted_neg_binomial_comparison_c4 = list(dist = "shifted_neg_binomial",
                                            truth_dist = "shifted_poisson", rate = c(5, 10),
                                            shift = c(10, 30), grid = 0:200, modes = 2L,
                                            title = "C4: Poisson control",
                                            fit = list(count = 500L, K = 6L,
                                                       nb_iter = 5000L, burnin = 2500L)),
  shifted_neg_binomial_comparison_c6 = list(dist = "shifted_neg_binomial",
                                            truth_dist = "zero_inflated_neg_binomial",
                                            zero_probability = 0.35, mean = 8, size = 2,
                                            grid = 0:200, modes = 2L,
                                            title = "C6: Excess zeros",
                                            fit = list(count = 500L, K = 6L,
                                                       nb_iter = 5000L, burnin = 2500L))
)

simulation_density <- function(experiment, x) {
  distribution <- if (is.null(experiment$truth_dist)) experiment$dist else experiment$truth_dist
  if (distribution == "zero_inflated_neg_binomial") {
    return(experiment$zero_probability * (x == 0) +
             (1 - experiment$zero_probability) * dnbinom(x, size = experiment$size,
                                                        mu = experiment$mean))
  }
  components <- vapply(1:2, function(component) {
    switch(distribution,
      normal = dnorm(x, experiment$location[component], experiment$scale),
      skew_normal = sn::dsn(x, experiment$location[component], experiment$scale,
                            experiment$shape[component]),
      poisson = dpois(x, experiment$rate[component]),
      shifted_poisson = dpois(x - experiment$shift[component], experiment$rate[component]),
      neg_binomial = dnbinom(x, size = experiment$size, mu = experiment$mean[component]),
      shifted_neg_binomial = dnbinom(x - experiment$shift[component],
                                     size = experiment$size, mu = experiment$mean[component])
    )
  }, numeric(length(x)))
  rowMeans(components)
}

simulation_sample <- function(experiment, count) {
  distribution <- if (is.null(experiment$truth_dist)) experiment$dist else experiment$truth_dist
  if (distribution == "zero_inflated_neg_binomial") {
    zeros <- runif(count) < experiment$zero_probability
    observations <- rnbinom(count, size = experiment$size, mu = experiment$mean)
    observations[zeros] <- 0L
    return(observations)
  }
  component <- sample.int(2, count, replace = TRUE)
  switch(distribution,
    normal = rnorm(count, experiment$location[component], experiment$scale),
    skew_normal = vapply(component, function(selected) {
      as.numeric(sn::rsn(1, experiment$location[selected], experiment$scale,
                         experiment$shape[selected]))
    }, numeric(1)),
    poisson = rpois(count, experiment$rate[component]),
    shifted_poisson = rpois(count, experiment$rate[component]) + experiment$shift[component],
    neg_binomial = rnbinom(count, size = experiment$size, mu = experiment$mean[component]),
    shifted_neg_binomial = rnbinom(count, size = experiment$size,
                                   mu = experiment$mean[component]) + experiment$shift[component]
  )
}

simulation_modes <- function(grid, density) {
  runs <- rle(density)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  peaks <- which(runs$values > 0 &
                   runs$values > c(-Inf, head(runs$values, -1)) &
                   runs$values > c(tail(runs$values, -1), -Inf))
  (grid[starts[peaks]] + grid[ends[peaks]]) / 2
}

mixture_density_from_draw <- function(draw, fit, x) {
  components <- lapply(seq_len(fit$K), function(component) {
    parameters <- vapply(fit$pars_names, function(name) {
      draw[[paste0(name, component)]]
    }, numeric(1))
    parameters[["eta"]] * fit$pdf_func(x, parameters)
  })
  Reduce(`+`, components)
}

test_that("simulation references distinguish components, modes, plateaus, and tails", {
  expect_equal(simulation_modes(0:6, c(0, 1, 1, 0, 2, 0, 0)), c(1.5, 4))
  expect_equal(simulation_modes(0:3, c(1, 0, 0, 0)), 0)
  for (experiment in simulation_experiments) {
    density <- simulation_density(experiment, experiment$grid)
    step <- diff(experiment$grid)[1]
    expect_equal(sum(density) * step, 1, tolerance = 1e-7)
    expect_length(simulation_modes(experiment$grid, density), experiment$modes)
    set.seed(42)
    observations <- simulation_sample(experiment, 10000)
    expect_true(is.vector(observations))
    expect_length(observations, 10000L)
    expect_true(all(is.finite(observations)))
    for (probability in c(0.1, 0.5, 0.9)) {
      index <- which(cumsum(density) * step >= probability)[1]
      cutoff <- experiment$grid[index]
      expected_cdf <- sum(density[experiment$grid <= cutoff]) * step
      expect_lt(abs(mean(observations <= cutoff) - expected_cdf), 0.025)
    }
  }
})

plot_directory <- Sys.getenv("BAYESMULTIMODE_SIMULATION_PLOT_DIR")
plot_results <- new.env(parent = emptyenv())
if (nzchar(plot_directory)) {
  dir.create(plot_directory, recursive = TRUE, showWarnings = FALSE)
  plot_script <- testthat::test_path("simulation", "plot-estimation.R")
  if (!file.exists(plot_script)) plot_script <- testthat::test_path("plot-estimation.R")
  source(plot_script, local = TRUE)
}
selected_method <- Sys.getenv("BAYESMULTIMODE_SIMULATION_METHOD")
experiment_names <- names(simulation_experiments)
if (nzchar(selected_method)) {
  experiment_names <- experiment_names[vapply(simulation_experiments, function(experiment) {
    experiment$dist == selected_method
  }, logical(1))]
  if (!length(experiment_names)) stop("Unknown simulation method: ", selected_method)
}

for (experiment_name in experiment_names) {
  experiment <- simulation_experiments[[experiment_name]]
  experiment_index <- match(experiment_name, names(simulation_experiments))
  fit_settings <- utils::modifyList(
    list(count = 600L, K = 10L, nb_iter = 2000L, burnin = 1000L),
    if (is.null(experiment$fit)) list() else experiment$fit
  )
  chain_summaries <- new.env(parent = emptyenv())
  chain_summaries$draws <- vector("list", 4)
  chain_summaries$densities <- vector("list", 4)
  chain_summaries$errors <- rep("", 4)
  run_convergence <- Sys.getenv("BAYESMULTIMODE_RUN_CONVERGENCE") == "true"

  for (chain in seq_along(chain_summaries$draws)) {
    test_that(paste(experiment_name, "chain", chain, "recovers the true distribution"), {
      skip_if(Sys.getenv("BAYESMULTIMODE_RUN_SIMULATION") != "true",
              "Set BAYESMULTIMODE_RUN_SIMULATION=true to run recovery checks.")
      grid <- experiment$grid
      step <- diff(grid)[1]
      truth <- simulation_density(experiment, grid)
      true_modes <- simulation_modes(grid, truth)
      set.seed(202600 + experiment_index)
      observations <- simulation_sample(experiment, fit_settings$count)
      chain_summaries$observations <- observations
      set.seed(10000 * experiment_index + chain)
      if (nzchar(plot_directory)) cat("\nFitting", experiment_name, "chain", chain, "\n")
      fit <- tryCatch(
        bayes_fit(observations, K = fit_settings$K, dist = experiment$dist,
                  nb_iter = fit_settings$nb_iter, burnin = fit_settings$burnin, print = FALSE),
        error = function(error) {
          chain_summaries$errors[chain] <- conditionMessage(error)
          stop(error)
        }
      )
      expect_s3_class(fit, "bayes_mixture")
      expect_equal(nrow(fit$mcmc_all), fit_settings$nb_iter)
      expect_equal(nrow(fit$mcmc), fit_settings$nb_iter - fit_settings$burnin)
      expect_true(all(is.finite(fit$mcmc)))
      expect_true(all(is.finite(tail(fit$loglik, nrow(fit$mcmc)))))
      weights <- fit$mcmc[, paste0("eta", seq_len(fit$K)), drop = FALSE]
      expect_equal(fit$K, fit_settings$K)
      expect_true(all(weights >= 0))
      expect_equal(rowSums(weights), rep(1, nrow(weights)), tolerance = 1e-10)

      density_draws <- vapply(seq_len(nrow(fit$mcmc)), function(draw) {
        mixture_density_from_draw(fit$mcmc[draw, ], fit, grid)
      }, numeric(length(grid)))
      mass <- colSums(density_draws) * step
      expect_true(all(is.finite(density_draws)))
      expect_true(all(abs(mass - 1) < 0.01))
      posterior_density <- rowMeans(density_draws)
      chain_summaries$densities[[chain]] <- posterior_density
      total_variation <- (sum(abs(posterior_density - truth)) * step +
                            abs(sum(posterior_density - truth) * step)) / 2
      expect_lt(total_variation, 0.15)
      cdf_error <- max(abs(cumsum(posterior_density - truth) * step))
      expect_lt(cdf_error, 0.1)

      if (run_convergence) {
        quantile_index <- vapply(c(0.1, 0.5, 0.9), function(probability) {
          which(cumsum(truth) * step >= probability)[1]
        }, integer(1))
        mean_draws <- as.vector(crossprod(grid, density_draws)) * step
        variance_draws <- as.vector(crossprod(grid^2, density_draws)) * step - mean_draws^2
        cdf_draws <- vapply(grid[quantile_index], function(probe) {
          colSums(density_draws[grid <= probe, , drop = FALSE]) * step
        }, numeric(nrow(fit$mcmc)))
        colnames(cdf_draws) <- c("cdf10", "cdf50", "cdf90")
        chain_summaries$draws[[chain]] <- cbind(
          loglik = tail(fit$loglik, nrow(fit$mcmc)),
          mean = mean_draws, sd = sqrt(variance_draws), cdf_draws
        )
      }

      selected <- unique(round(seq(1, nrow(fit$mcmc), length.out = 100)))
      recalculated <- vapply(selected, function(draw) {
        sum(log(mixture_density_from_draw(fit$mcmc[draw, ], fit, observations)))
      }, numeric(1))
      expect_equal(recalculated, tail(fit$loglik, nrow(fit$mcmc))[selected], tolerance = 1e-8)
      mode_input <- fit
      mode_input$mcmc <- fit$mcmc[selected, , drop = FALSE]
      modes <- bayes_mode(mode_input, rd = 2, tol_x = 0.1)
      counts <- rowSums(!is.na(modes$modes))
      expect_gt(mean(counts == length(true_modes)), 0.7)
      location_tolerance <- if (step == 1) 2 else 0.5
      recovered <- apply(modes$modes, 1, function(locations) {
        locations <- sort(locations[is.finite(locations)])
        length(locations) == length(true_modes) &&
          all(abs(locations - true_modes) <= location_tolerance)
      })
      expect_gt(mean(recovered), 0.6)
    })
  }

  if (nzchar(plot_directory) && Sys.getenv("BAYESMULTIMODE_RUN_SIMULATION") == "true") {
    plot_results[[experiment_name]] <- list(
      experiment = experiment, observations = chain_summaries$observations,
      densities = chain_summaries$densities, errors = chain_summaries$errors,
      fit = fit_settings
    )
  }

  test_that(paste(experiment_name, "has mixture-level R-hat and ESS diagnostics"), {
    skip_if(!run_convergence, "Mixture-level convergence checks are disabled.")
    skip_if(Sys.getenv("BAYESMULTIMODE_RUN_SIMULATION") != "true",
            "Enable recovery simulations to generate the diagnostic chains.")
    expect_false(any(vapply(chain_summaries$draws, is.null, logical(1))))
    quantities <- if (any(vapply(chain_summaries$draws, is.null, logical(1)))) {
      character()
    } else colnames(chain_summaries$draws[[1]])
    for (quantity in quantities) {
      chains <- vapply(chain_summaries$draws, function(summary) summary[, quantity],
                       numeric(nrow(chain_summaries$draws[[1]])))
      diagnostics <- c(rhat = posterior::rhat(chains),
                       bulk_ess = posterior::ess_bulk(chains),
                       tail_ess = posterior::ess_tail(chains))
      expect_true(all(is.finite(diagnostics)), info = quantity)
      expect_true(diagnostics[["rhat"]] < 1.05,
                  info = paste(quantity, "R-hat =", round(diagnostics[["rhat"]], 3)))
      expect_true(diagnostics[["bulk_ess"]] > 400,
                  info = paste(quantity, "bulk ESS =", round(diagnostics[["bulk_ess"]])))
      expect_true(diagnostics[["tail_ess"]] > 400,
                  info = paste(quantity, "tail ESS =", round(diagnostics[["tail_ess"]])))
    }
  })
}

if (nzchar(plot_directory) && Sys.getenv("BAYESMULTIMODE_RUN_SIMULATION") == "true") {
  methods <- unique(vapply(simulation_experiments[experiment_names], function(experiment) {
    experiment$dist
  }, character(1)))
  for (method in methods) {
    method_names <- experiment_names[vapply(simulation_experiments[experiment_names], function(experiment) {
      experiment$dist == method
    }, logical(1))]
    method_results <- mget(method_names, envir = plot_results)
    plot_simulation_method(method_results, file.path(plot_directory, paste0(method, ".png")))
  }
}
