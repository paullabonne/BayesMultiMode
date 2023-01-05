#' Bayesian estimation of mixture distributions
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numeric vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
bayes_estimation <- function(data,
                             K,
                             dist,
                             nb_iter = 2000,
                             burnin = nb_iter/2,
                             chains = 4,
                             cores = 4,
                             a0 = 10,
                             A0 = 10*K,
                             b0 = median(data),
                             B0 = (max(data) - min(data))^2,
                             c0 = 2.5,#2.5,
                             e0 = 0,
                             g0 = 0.5,
                             G0 = 100*2.5/0.5/(max(data) - min(data))^2,#0.5/(sd(y)^2/2),
                             #skew prior
                             h0 = 0,
                             H0 = 10,
                             #studen t prior
                             n0 = 2,
                             N0 = 0.1,
                             refresh = 1e3
) {

  mixture_data <- list(K = K,
                       N = length(data),
                       y = data,
                       a0 = a0,
                       A0 = A0,
                       b0 = b0,
                       B0 = B0,
                       c0 = c0,
                       e0 = e0,
                       g0 = g0,
                       G0 = G0)

  if (dist %in% c("skew_t", "skew_normal")) {
    mixture_data$h0 = h0
    mixture_data$H0 = H0
  }

  if (dist %in% c("student", "skew_t")) {
    mixture_data$n0 = n0
    mixture_data$N0 = N0
  }

  # if (dist == "skew_normal") {
  #   mix_model = stanmodels$skew_normal_mixture
  # }
  # if (dist == "gaussian") {
  #   mix_model = stanmodels$gaussian_mixture
  # }
  # if (dist == "student") {
  #   mix_model = stanmodels$student_mixture
  # }
  # if (dist == "skew_r") {
  #   mix_model = stanmodels$skew_t_mixture
  # }

  fit <- rstan::sampling(stanmodels[[paste0(dist, "_mixture")]],  # Stan program
                         data = mixture_data,    # named list of data
                         chains = chains,             # number of Markov chains
                         warmup = burnin,        # number of warmup iterations per chain
                         iter = nb_iter,         # total number of iterations per chain
                         cores = cores,              # number of cores (could use one per chain)
                         refresh = refresh           # no progress shown
  )

  return(fit)
}
