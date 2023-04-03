#' Bayesian estimation of mixture distributions
#' 
#' Gibbs samplers for sparse finite mixture Markov chain Monte Carlo (SFM MCMC) estimation.
#' 
#' @param data Vector of observations
#' @param K Maximum number of mixture components
#' @param dist String indicating the distribution of the mixture components
#' Currently supports "normal", "skew_normal", "poisson" and "shifted_poisson"
#' @param priors List of priors; default is an empty list which implies the following priors :\cr
#' a0 = 1,\cr A0 = 200,\cr b0 = median(y),\cr B0 = (max(y) - min(y))^2 (normal),\cr
#' D_xi = 1,\cr D_psi =1, (skew normal: B0 = diag(D_xi,D_psi)), \cr c0 = 2.5,\cr
#' l0 = 1.1 (poisson),\cr l0 = 5 (shifted poisson),\cr L0 = 1.1/median(y),\cr L0 = l0 - 1 (shifted poisson),\cr
#' g0 = 0.5,\cr G0 = 100*g0/c0/B0 (normal),\cr 
#' G0 = g0/(0.5*var(y)) (skew normal)
#' @param nb_iter Number of MCMC iterations; default is 2000
#' @param burnin Number of MCMC iterations used as burnin; default is nb_iter/2
#' @param printing Showing MCMC progression ?
#' 
#' @return A list of class `BayesMixture` containing
#' \itemize{
#'  \item{data}{ - Same as argument}
#'  \item{dist_type}{ - Type of the distribution (continuous or discrete)}
#'  \item{pars_names}{ - Names of the mixture components' parameters}
#'  \item{mcmc}{ - Matrix of MCMC draws where the rows corresponding to burnin have been discarded}
#'  \item{mcmc_all}{ - Original matrix of MCMC draws}
#' }
#' 
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr\cr
#' \insertRef{fruhwirth-schnatter_bayesian_2010}{BayesMultiMode}\cr\cr
#' \insertRef{SFS:Mal:2019}{BayesMultiMode}\cr\cr
#' \insertRef{viallefont2002bayesian}{BayesMultiMode}
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' 
#' @examples
#' # Example with galaxy data ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = galaxy
#'
#' # estimation
#' bayesmix = bayes_estimation(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'                            
#' # plot estimated mixture
#' # plot(bayesmix, max_size = 200)
#' 
#' # Example with DNA data =====================================================
#' \donttest{
#' set.seed(123) 
#' 
#' # retrieve DNA data
#' y = d4z4
#'
#' # estimation
#' bayesmix = bayes_estimation(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "shifted_poisson",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'                            
#' # plot estimated mixture
#' # plot(bayesmix, max_size = 200)
#' }
#' 
#' @export
bayes_estimation <- function(data,
                             K,
                             dist,
                             priors = list(),
                             nb_iter = 2000,
                             burnin = nb_iter/2,
                             printing = TRUE) {
  
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(!any(is.na(data)) & !any(is.infinite(data)),
              msg = "y should not include missing or infinite values")
  assert_that(dist %in% c("normal", "skew_normal", "poisson", "shifted_poisson") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, skew_normal, poisson, shifted_poisson or 'NA'")
  assert_that(is.scalar(nb_iter) & nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin) & burnin > 0 & burnin < nb_iter,
              msg = "nb_iter should be a positive integer lower than burnin")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  assert_that(is.logical(printing), msg = "printing should be either TRUE or FALSE")

  # rounding parameters that should be integers
  K = round(K)
  nb_iter = round(nb_iter)
  burnin = round(burnin)
  
  if (dist %in% c("normal")) {
    priors_labels = c("a0", "A0", "e0", "b0", "B0", "c0", "g0", "G0")
    
    mcmc = gibbs_SFM_normal(y = data,
                            K = K,
                            nb_iter = nb_iter,
                            priors = priors[priors_labels],
                            printing = printing)
    pars_names = c("eta", "mu", "sigma")
    dist_type = "continuous"
    
  } else if (dist == "skew_normal") {
    priors_labels = c("a0", "A0", "e0", "b0", "c0", "C0", "g0", "G0", "D_xi", "D_psi")
    
    mcmc <- gibbs_SFM_skew_n(y = data,
                             K = K,
                             nb_iter = nb_iter,
                             priors = priors[priors_labels],
                             printing = printing)
    pars_names = c("eta", "xi", "omega", "alpha")
    dist_type = "continuous"
    
  } else if (dist == "poisson") {
    priors_labels = c("a0", "A0", "e0", "l0", "L0")
    
    mcmc <- gibbs_SFM_poisson(y = data,
                              K = K,
                              nb_iter = nb_iter,
                              priors = priors[priors_labels],
                              printing = printing)
    pars_names = c("eta", "lambda")
    dist_type = "discrete"
    
  } else if (dist == "shifted_poisson") {
    priors_labels = c("a0", "A0", "e0", "l0", "L0")
    
    mcmc <- gibbs_SFM_sp(y = data,
                         K = K,
                         nb_iter = nb_iter,
                         priors = priors[priors_labels],
                         printing = printing)
    pars_names = c("eta", "kappa", "lambda")
    dist_type = "discrete"
    
  } else {
    stop("mixture distribution not supported")
  }
  
  BayesMixture = new_BayesMixture(mcmc = mcmc,
                                  data = data, K = K,
                                  burnin = burnin, dist = dist,
                                  dist_type = dist_type,
                                  pars_names = pars_names)
  
  return(BayesMixture)
}