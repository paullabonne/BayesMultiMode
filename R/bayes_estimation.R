#' Bayesian estimation of mixture distributions
#'
#' @param data Numeric vector of input values.
#' @param K Integer indicating the maximum number of mixture components.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "student" and "skew_normal".
#' @param fit Numeric vector of output values.
#' @param e0 Numeric positive value for the dirichlet prior parameters. If 0 so an hyperprior is used instead. Default is 0.
#' @param a0 Numeric value for the means of the gamma hyperprior used the dirichlet prior. Default is 10.
#' @param A0 Numeric value for the variance of the gamma hyperprior used the dirichlet prior. Default is 10*K.
#' @param b0 Numeric value for the means of the mean priors. Default is mean(data).
#' @param B0 Numeric value for the variance of the mean priors. Default is R^2 where R = (max(data) - min(data)).
#' @param c0 Numeric value for variance prior. Default is 2.5.
#' @param g0 Numeric value for variance inverse gamma hyperprior. Default is 0.5.
#' @param G0 Numeric value for variance inverse gamma hyperprior. Default is 100*2.5/0.5/R^2.
#' @param h0 Numeric value for the means of the skew parameters priors. Default is 0.
#' @param H0 Numeric value for the variance of the skew parameters priors. Default is 10.
#' @param n0 Numeric value for the means of the degree of freedom gamma priors. Default is 2.
#' @param N0 Numeric value for the variance of the degree of freedom gamma priors. Default is 0.1.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' 
#' @importFrom rstan sampling
#' 
#' @return An object of class `BayesMixture`.
#' 
#' @export
#'
#'
bayes_estimation <- function(data,
                             K,
                             dist,
                             a0 = 10,
                             A0 = 10*K,
                             b0 = mean(data),
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
                             nb_iter = 2000,
                             burnin = nb_iter/2,
                             chains = 4,
                             cores = 4,
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
  
  if (dist %in% c("normal", "student", "skew_normal")) {
    fit <- sampling(stanmodels[[paste0(dist, "_mixture")]],  # Stan program
                    data = mixture_data,    # named list of data
                    chains = chains,             # number of Markov chains
                    warmup = burnin,        # number of warmup iterations per chain
                    iter = nb_iter,         # total number of iterations per chain
                    cores = cores,              # number of cores (could use one per chain)
                    refresh = refresh           # no progress shown
    ) 
    
    dist_type = "continuous"
    
  } else if (dist %in% c("shifted_poisson")) {
    fit <- shift_pois_mcmc(y = data, K, nb_iter, burnin)
    
    dist_type = "discrete"
    
  } else {
    stop("mixture distribution not supported")
  }
  
  BayesMixture = new_BayesMixture(fit, data, dist, dist_type = dist_type)

  return(BayesMixture)
}
