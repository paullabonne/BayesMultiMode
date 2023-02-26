#' Bayesian estimation of mixture distributions
#'
#' @param data Numeric vector of input values.
#' @param K Integer indicating the maximum number of mixture components.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "student" and "skew_normal".
#' @param e0 Dirichlet prior parameter. Must be positive or zero. If set to 0 then an hyperprior is used instead. Default is 0.
#' @param a0 Mean of the gamma hyperprior used for the dirichlet prior. Default is 10.
#' @param A0 Variance of the gamma hyperprior used for the dirichlet prior. Default is 10*K.
#' @param b0 Mean of the mean priors. Default is mean(data).
#' @param B0 Variance of the mean priors. Default is R^2 where R = (max(data) - min(data)).
#' @param c0 Variance prior. Default is 2.5.
#' @param g0 Variance inverse gamma hyperprior. Default is 0.5.
#' @param G0 Variance inverse gamma hyperprior. Default is 100*2.5/0.5/R^2.
#' @param h0 Mean of the skew parameters priors. Default is 0.
#' @param H0 Variance of the skew parameters priors. Default is 10.
#' @param n0 Mean of the degree of freedom gamma priors. Default is 2.
#' @param N0 Variance of the degree of freedom gamma priors. Default is 2.
#' @param l0 ... Prior for lambda 
#' @param L0 ... Prior for lambda 
#' @param e0_kappa Dirichlet prior parameter for kappa (shifted poisson). Must be positive or zero. If set to 0 then an hyperprior is used instead.
#' Default is 1e-20
#' @param d0 ...
#' @param D0 ...
#' @param nb_iter Number of MCMC iterations. Default is 2000.
#' @param burnin Number of MCMC iterations used as burnin.
#' @param chains Number of chains.
#' @param cores Number of cores for parallel computation.
#' @param refresh Show intemediate results. Default is nb_iter/10 so intermediate results are shown every 1000 iterations.
#' @param ... Other arguments passed to `rstan::sampling`.
#' 
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
                             c0 = 2.5,
                             e0 = 0,
                             g0 = 0.5,
                             G0 = 100*g0/c0/B0,
                             #skew prior
                             h0 = 0,
                             H0 = 10,
                             #studen t prior
                             n0 = 2,
                             N0 = 0.1,
                             #poisson prior
                             l0 = 1.1,
                             L0 = 1.1/mean(data),
                             #shifted poisson prior
                             e0_kappa = 1e-20,
                             d0 = 10,
                             D0 = 10*(max(data) - min(data)),
                             # mcmc stan pars
                             nb_iter = 2000,
                             burnin = nb_iter/2,
                             chains = 4,
                             cores = 4,
                             refresh = nb_iter/10,
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
  assert_that(dist %in% c("normal", "student", "skew_t", "shifted_poisson_bis",
                          "skew_normal", "poisson", "shifted_poisson") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, student,
              skew_normal, skew_t, poisson, shifted_poisson, shifted_poisson_bis or NA")
  assert_that(is.scalar(nb_iter) & nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin) & burnin > 0 & burnin < nb_iter,
              msg = "nb_iter should be a positive integer lower than burnin")
  assert_that(is.scalar(chains) & chains > 0, msg = "chains should be a positive integer")
  assert_that(is.scalar(cores) & cores > 0, msg = "cores should be a positive integer")
  assert_that(is.scalar(refresh) & refresh > 0, msg = "refresh should be a positive integer")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  
  assert_that(is.scalar(A0) & A0 > 0, msg = "A0 should be positive")
  assert_that(is.scalar(B0) & B0 > 0, msg = "B0 should be a positive integer")
  assert_that(is.scalar(H0) & H0 > 0, msg = "H0 should be a positive integer")
  assert_that(is.scalar(N0) & N0 > 0, msg = "N0 should be a positive integer")
  assert_that(is.scalar(L0) & L0 > 0, msg = "L0 should be a positive integer")
  assert_that(is.scalar(D0) & D0 > 0, msg = "D0 should be a positive integer")
  
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
                       G0 = G0,
                       h0 = h0,
                       H0 = H0,
                       n0 = n0,
                       N0 = N0,
                       l0 = l0,
                       L0 = L0,
                       e0_kappa = e0_kappa,
                       d0 = d0,
                       D0 = D0)
  
  if (dist == "skew_normal") {
    pars_names = c("theta", "mu", "sigma", "xi")
  }
  
  if (dist == "student") {
    pars_names = c("theta", "mu", "sigma", "nu")
  }
  
  if (dist == "skew_t") {
    pars_names = c("theta", "mu", "sigma", "xi", "nu")
  }
  
  if (dist %in% c("normal")) {
    fit = gibbs_SFM_normal(y = data,
                           K,
                           nb_iter)
    pars_names = c("theta", "mu", "sigma")
    dist_type = "continuous"
  } else if (dist == "shifted_poisson") {
    fit <- gibbs_SFM_sp(y = data, K, nb_iter)
    pars_names = c("theta", "kappa", "lambda")
    dist_type = "discrete"
    
  } else {
    stop("mixture distribution not supported")
  }
  
  
  attr(fit, "K") = K
  attr(fit, "warmup") = burnin
  
  BayesMixture = new_BayesMixture(fit, data, dist, dist_type = dist_type, pars_names = pars_names)
  
  return(BayesMixture)
}
