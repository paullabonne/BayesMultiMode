tele_log_prior_K <- function(total, prior) {
  if (prior == "uniform") return(rep(0, length(total)))
  failures <- total - 1
  lbeta(5, 3 + failures) - lbeta(4, 3)
}

tele_log_partition <- function(total, counts, concentration, model) {
  shape <- if (model == "dynamic") concentration / total else concentration
  mass <- if (model == "dynamic") concentration else total * concentration
  result <- lgamma(mass) - lgamma(mass + sum(counts))
  for (size in counts) result <- result + lgamma(shape + size) - lgamma(shape)
  result
}

tele_log_count <- function(total, counts, concentration, model, prior) {
  base <- tele_log_prior_K(total, prior) + lfactorial(total) - lfactorial(total - length(counts))
  if (model == "static") {
    mass <- total * concentration
    return(base + lbeta(mass, sum(counts)))
  }
  base + tele_log_partition(total, counts, concentration, model)
}

tele_log_concentration <- function(log_value, total, counts, model) {
  value <- exp(log_value)
  if (!is.finite(value) || value == 0) return(-Inf)
  log_prior <- if (model == "dynamic") stats::df(value, 6, 3, log = TRUE) else
    stats::dgamma(value, shape = 1, rate = 20, log = TRUE)
  tele_log_partition(total, counts, value, model) + log_prior + log_value
}

tele_normal_family <- function(data, priors) {
  empty <- function(size, rate) {
    cbind(mu = rnorm(size, priors$b0, sqrt(priors$B0)),
          sigma = sqrt(1 / rgamma(size, priors$c0, rate)))
  }
  list(fields = c("mu", "sigma"), initial_shared = priors$g0,
       initialize = function(size, rate) {
         pars <- empty(size, rate)
         pars[, "mu"] <- data[sample.int(length(data), size, replace = TRUE)]
         pars
       },
       log_density = function(values, pars) dnorm(values, pars["mu"], pars["sigma"], log = TRUE),
       occupied = function(values, pars, rate, iteration) {
         draw <- draw_normal_component(values, pars["mu"], rate, priors)
         c(draw[1], sqrt(draw[2]))
       },
       shared = function(pars) draw_normal_rate(pars[, "sigma"]^2, priors),
       empty = empty)
}

tele_skew_normal_family <- function(data, priors) {
  prior_precision <- diag(1 / c(priors$D_xi, priors$D_psi))
  prior_score <- prior_precision %*% c(priors$b0, 0)
  parameters <- function(xi, psi, variance) {
    cbind(xi = xi, psi = psi, variance = variance,
          omega = sqrt(variance + psi^2), alpha = psi / sqrt(variance))
  }
  empty <- function(size, rate) {
    parameters(rnorm(size, priors$b0, sqrt(priors$D_xi)), rnorm(size, 0, sqrt(priors$D_psi)),
               1 / rgamma(size, priors$c0, rate = rate))
  }
  list(fields = c("xi", "omega", "alpha"), initial_shared = priors$C0,
       initialize = function(size, rate) {
         pars <- empty(size, rate)
         pars[, "xi"] <- data[sample.int(length(data), size, replace = TRUE)]
         pars
       },
       log_density = function(values, pars) {
         dsn(values, pars[["xi"]], pars[["omega"]], pars[["alpha"]], log = TRUE)
       },
       occupied = function(values, pars, rate, iteration) {
         design <- skew_normal_design(values, pars[["xi"]], pars[["psi"]], pars[["variance"]])
         residuals <- values - design %*% pars[c("xi", "psi")]
         variance <- 1 / rgamma(1, priors$c0 + length(values) / 2, rate + sum(residuals^2) / 2)
         regression <- skew_normal_regression(design, values, variance, prior_precision, prior_score)
         coefficients <- rmvnorm(1, regression$mean, regression$covariance)
         parameters(coefficients[1], coefficients[2], variance)[1, ]
       },
       shared = function(pars) draw_normal_rate(pars[, "variance"], priors),
       empty = empty)
}

tele_poisson_family <- function(data, priors, shifted) {
  bound <- max(data)
  initialize <- function(size, shared) {
    cbind(kappa = 0L, lambda = rgamma(size, priors$l0, rate = priors$L0))
  }
  list(fields = c(if (shifted) "kappa", "lambda"), initialize = initialize,
       log_density = function(values, pars) {
         dpois_unique(values - pars[["kappa"]], pars[["lambda"]], log = TRUE)
       },
       occupied = function(values, pars, shared, iteration) {
         shift <- if (shifted) draw_kap_lambda(values, pars[["kappa"]], pars[["lambda"]],
                                               priors$l0, priors$L0)$kappa else 0L
         c(kappa = shift, lambda = rgamma(1, priors$l0 + sum(values - shift),
                                          rate = priors$L0 + length(values)))
       },
       empty = function(size, shared) {
         pars <- initialize(size, shared)
         if (shifted) pars[, "kappa"] <- sample.int(bound + 1, size, replace = TRUE) - 1L
         pars
       })
}

tele_nb_family <- function(data, priors, burnin, shifted) {
  bound <- max(data)
  if (bound <= 0) stop("The negative binomial mean prior requires a positive full-data maximum.")
  if (!is.numeric(priors$sig) || length(priors$sig) != 1L || priors$sig <= 0)
    stop("sig must be a positive proposal scale.")
  empty <- function(size, shared) {
    draws <- replicate(size, rneg_binom(numeric(), priors$sig, 0, 0, bound))
    cbind(kappa = if (shifted) sample.int(bound + 1, size, replace = TRUE) - 1L else 0L,
          mu = draws[1, ], rtilde = draws[2, ], scale = priors$sig, accepted = 0, proposed = 0)
  }
  list(fields = c(if (shifted) "kappa", "mu", "rtilde"),
       initialize = function(size, shared) {
         pars <- empty(size, shared)
         pars[, "kappa"] <- 0L
         pars
       },
       log_density = function(values, pars) {
         dnbinom_unique(values - pars["kappa"], size = (1 - pars["rtilde"]) / pars["rtilde"],
                        mu = pars["mu"], log = TRUE)
       },
       occupied = function(values, pars, shared, iteration) {
         if (shifted) {
           shift <- draw_kap_mu(values, pars[["kappa"]], pars[["mu"]], pars[["rtilde"]], bound)
           pars[c("kappa", "mu")] <- c(shift$kappa, shift$mu)
         }
         draw <- rneg_binom(values - pars[["kappa"]], pars[["scale"]], pars[["mu"]], pars[["rtilde"]], bound)
         pars[c("mu", "rtilde")] <- draw[1:2]
         if (iteration <= burnin) {
           pars["accepted"] <- pars["accepted"] + draw[3]
           pars["proposed"] <- pars["proposed"] + 1
           if (iteration %% 50 == 0) {
             pars["scale"] <- adapt_rtilde_proposals(pars["scale"], pars["accepted"], pars["proposed"], iteration / 50)
             pars[c("accepted", "proposed")] <- 0
           }
         }
         pars
       },
       empty = empty)
}

tele_density <- function(data, pars, weights, family) {
  log_terms <- matrix(0, length(data), nrow(pars))
  for (component in seq_len(nrow(pars))) {
    log_terms[, component] <- log(weights[component]) + family$log_density(data, pars[component, ])
  }
  maximum <- log_terms[, 1]
  for (component in seq_len(ncol(log_terms))[-1]) maximum <- pmax(maximum, log_terms[, component])
  if (any(!is.finite(maximum))) stop("No finite mixture density for an observation.")
  relative <- exp(log_terms - maximum)
  list(relative = relative, loglik = sum(maximum + log(rowSums(relative))))
}

tele_occupied <- function(allocations, pars) {
  counts <- colSums(allocations)
  occupied <- which(counts > 0)
  list(counts = counts[occupied], pars = pars[occupied, , drop = FALSE],
       allocations = allocations[, occupied, drop = FALSE])
}

gibbs_telescoping <- function(data, cap, initial, iterations, family, model,
                              prior, fixed, print) {
  fields <- c("eta", family$fields)
  columns <- unlist(lapply(fields, paste0, seq_len(cap)), use.names = FALSE)
  draws <- matrix(NA_real_, iterations, length(columns), dimnames = list(NULL, columns))
  state <- matrix(NA_real_, iterations, 3, dimnames = list(NULL, c("K", "Kplus", "concentration")))
  loglik <- numeric(iterations)
  concentration <- if (is.null(fixed)) 1 else fixed
  # Log-scale proposal SDs from telescope's F_6_3 and G_1_20 specifications.
  proposal_sd <- if (model == "dynamic") 2.5 else 1.5
  shared <- family$initial_shared
  pars <- family$initialize(initial, shared)
  density <- tele_density(data, pars, rep(1 / initial, initial), family)
  for (iteration in seq_len(iterations)) {
    partition <- tele_occupied(draw_allocations(density$relative), pars)
    counts <- partition$counts
    occupied <- length(counts)
    pars <- partition$pars
    for (component in seq_len(occupied)) {
      values <- data[partition$allocations[, component] == 1L]
      pars[component, ] <- family$occupied(values, pars[component, ], shared, iteration)
    }
    if (!is.null(family$shared)) shared <- family$shared(pars)
    support <- seq.int(occupied, cap)
    log_mass <- tele_log_count(support, counts, concentration, model, prior)
    total <- support[sample.int(length(support), 1, prob = exp(log_mass - max(log_mass)))]
    if (is.null(fixed)) {
      current <- log(concentration)
      proposal <- rnorm(1, current, proposal_sd)
      ratio <- tele_log_concentration(proposal, total, counts, model) -
        tele_log_concentration(current, total, counts, model)
      if (log(runif(1)) < ratio) concentration <- exp(proposal)
    }
    shape <- if (model == "dynamic") concentration / total else concentration
    if (total > occupied) pars <- rbind(pars, family$empty(total - occupied, shared))
    weights <- as.numeric(rdirichlet(1, c(counts, numeric(total - occupied)) + shape))
    density <- tele_density(data, pars, weights, family)
    saved <- cbind(eta = weights, pars[, family$fields, drop = FALSE])
    if (any(!is.finite(saved)) || abs(sum(weights) - 1) > 1e-10 ||
        !is.finite(density$loglik) || sum(counts) != length(data)) stop("Invalid telescoping draw.")
    slots <- outer(seq_len(total), (seq_along(fields) - 1L) * cap, `+`)
    draws[iteration, as.vector(slots)] <- as.vector(saved)
    state[iteration, ] <- c(total, occupied, concentration)
    loglik[iteration] <- density$loglik
    if (print && iteration %% max(1, round(iterations / 10)) == 0)
      cat(100 * iteration / iterations, "% draws finished\n")
  }
  list(draws = draws, state = state, loglik = loglik)
}