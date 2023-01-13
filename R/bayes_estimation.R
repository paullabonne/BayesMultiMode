#' Bayesian estimation of mixture distributions
#'
#' @param data Numeric vector of input values.
#' @param K Integer indicating the maximum number of mixture components.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "student" and "skew_normal".
#' @param e0 (numeric positive scalar) Dirichlet prior parameter. If 0 so an hyperprior is used instead. Default is 0.
#' @param a0 (numeric scalar) Mean of the gamma hyperprior used the dirichlet prior. Default is 10.
#' @param A0 (numeric scalar) Variance of the gamma hyperprior used the dirichlet prior. Default is 10*K.
#' @param b0 (numeric scalar) Mean of the mean priors. Default is mean(data).
#' @param B0 (numeric scalar) Numeric value for the variance of the mean priors. Default is R^2 where R = (max(data) - min(data)).
#' @param c0 (numeric scalar) Numeric value for variance prior. Default is 2.5.
#' @param g0 (numeric scalar) Numeric value for variance inverse gamma hyperprior. Default is 0.5.
#' @param G0 (numeric scalar) Numeric value for variance inverse gamma hyperprior. Default is 100*2.5/0.5/R^2.
#' @param h0 (numeric scalar) Mean of the skew parameters priors. Default is 0.
#' @param H0 (numeric scalar) Variance of the skew parameters priors. Default is 10.
#' @param n0 (numeric scalar) Mean of the degree of freedom gamma priors. Default is 2.
#' @param N0 (numeric scalar) Variance of the degree of freedom gamma priors. Default is 2.
#' @param nb_iter Number of MCMC iterations. Default is 2000.
#' @param burnin ...
#' @param chains ...
#' @param cores ...
#' @param refresh ...
#' @param ... Other arguments passed to `rstan::sampling`.
#' 
#' @importFrom rstan sampling
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar

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
                             refresh = 1e3,
                             ...
) {
  K = round(K)
  nb_iter = round(nb_iter)
  burnin = round(burnin)
  refresh = round(refresh)
  cores = round(cores)
  chains = round(chains)
  
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(dist %in% c("normal", "student", "skew_normal", "shifted_poisson") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, student, skew_normal, shifted_poisson or NA")
  assert_that(is.scalar(nb_iter) & nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin) & burnin > 0 & burnin < nb_iter,
              msg = "nb_iter should be a positive integer lower than burnin")
  assert_that(is.scalar(chains) & chains > 0, msg = "chains should be a positive integer")
  assert_that(is.scalar(cores) & cores > 0, msg = "cores should be a positive integer")
  assert_that(is.scalar(refresh) & refresh > 0, msg = "refresh should be a positive integer")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  

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
                    refresh = refresh,           # no progress shown
                    ...
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
