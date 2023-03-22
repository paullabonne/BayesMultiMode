#' Bayesian estimation of mixture distributions
#' 
#' Gibbs sampler for Spare Finite Mixture MCMC estimation of mixture distributions.
#' 
#' @param data Numeric vector of input values.
#' @param K Integer indicating the maximum number of mixture components.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "skew_normal", "poisson" and "shifted_poisson".
#' @param priors List of priors. Default in an empty list.
#' @param nb_iter Number of MCMC iterations. Default is 2000.
#' @param burnin Number of MCMC iterations used as burnin. Default is nb_iter/2.
#' @param printing Showing MCMC progression ?
#' 
#' @return An object of class `BayesMixture`.
#' 
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr\cr
#' \insertRef{SFS:Mal:2019}{BayesMultiMode}\cr
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' 
#' @export
#'
#'
bayes_estimation <- function(data,
                             K,
                             dist,
                             priors = list(),
                             nb_iter = 2000,
                             burnin = nb_iter/2,
                             printing = TRUE) {
  K = round(K)
  nb_iter = round(nb_iter)
  burnin = round(burnin)
  
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(dist %in% c("normal", "skew_normal", "poisson", "shifted_poisson") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, skew_normal, poisson, shifted_poisson or NA")
  assert_that(is.scalar(nb_iter) & nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin) & burnin > 0 & burnin < nb_iter,
              msg = "nb_iter should be a positive integer lower than burnin")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  
  if (dist %in% c("normal")) {
    priors_labels = c("a0", "A0", "e0", "b0", "B0", "c0", "g0", "G0")
    
    mcmc = gibbs_SFM_normal(y = data,
                            K = K,
                            nb_iter = nb_iter,
                            priors = priors[priors_labels],
                            printing = printing)
    pars_names = c(theta = "theta", mu = "mu", sigma = "sigma")
    dist_type = "continuous"
    
  } else if (dist == "skew_normal") {
    priors_labels = c("a0", "A0", "e0", "b0", "c0", "C0", "g0", "G0", "D_xi", "D_psi")
    
    mcmc <- gibbs_SFM_skew_n(y = data,
                             K = K,
                             nb_iter = nb_iter,
                             priors = priors[priors_labels],
                             printing = printing)
    pars_names = c(theta = "theta", xi = "xi", omega = "omega", alpha = "alpha")
    dist_type = "continuous"
    
  } else if (dist == "poisson") {
    priors_labels = c("a0", "A0", "e0", "l0", "L0")
    
    mcmc <- gibbs_SFM_poisson(y = data,
                              K = K,
                              nb_iter = nb_iter,
                              priors = priors[priors_labels],
                              printing = printing)
    pars_names = c(theta = "theta", lambda = "lambda")
    dist_type = "discrete"
    
  } else if (dist == "shifted_poisson") {
    priors_labels = c("a0", "A0", "e0", "l0", "L0")
    
    mcmc <- gibbs_SFM_sp(y = data,
                         K = K,
                         nb_iter = nb_iter,
                         priors = priors[priors_labels],
                         printing = printing)
    pars_names = c(theta = "theta", kappa = "kappa", lambda = "lambda")
    dist_type = "discrete"
    
  } else {
    stop("mixture distribution not supported")
  }
  
  BayesMixture = new_BayesMixture(mcmc, data, K, burnin, dist, dist_type = dist_type, pars_names = pars_names)
  
  return(BayesMixture)
}