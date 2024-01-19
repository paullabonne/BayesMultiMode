#' Creating a S3 object of class \code{BayesMixture}
#' 
#' Function for creating an object of class \code{BayesMixture} which can subsequently be used as argument in [bayes_mode()].
#' This function is useful for users who want to use the mode inference functions of the package with MCMC output generated using 
#' other software packages.
#' 
#' @param mcmc A matrix of MCMC draws.
#' @param data A vector containing the data used for estimating the model and generating the MCMC draws.
#' @param K Number of mixture components.
#' @param burnin Number of draws to discard as burnin.
#' @param dist Distribution family of the mixture components supported by
#' the package (e.g. "normal", "student", "skew_normal" or "shifted_poisson").
#' @param pdf_func Pdf or pmf of the mixture components;
#' this input is used only if dist_name is invalid or NULL.
#' @param dist_type Either "continuous" or "discrete".
#' @param loglik Vector showing the log likelihood at each MCMC draw.
#' 
#' @returns
#' A list of class \code{BayesMixture} containing:
#' \itemize{
#'  \item{data}{ - Same as argument.}
#'  \item{dist_type}{ - Same as argument.}
#'  \item{pars_names}{ - Names of distribution parameters.}
#'  \item{mcmc}{ - Matrix of MCMC draws where the rows corresponding to burnin have been discarded.}
#'  \item{mcmc_all}{ - Original matrix of MCMC draws.}
#' }
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom posterior subset_draws
#' @importFrom stringr str_extract
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_replace
#' @importFrom stringr str_locate
#' 
#' @examples
#' 
#' # Example with a Student t ================================================
#' mu = c(0.5,6)
#' sigma = c(1,2)
#' nu = c(5,5)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = mu, sigma = sigma, nu = nu)
#' pars_names = c("eta", "mu", "sigma", "nu")
#' dist_type = "continuous"
#'
#' data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
#'          sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
#'
#' fit = c(eta = p, mu = mu, sigma = sigma, nu = nu)
#' fit = rbind(fit, fit)
#' 
#' pdf_func = function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' BM = new_BayesMixture(fit, data, K = 2, burnin = 1, pdf_func = pdf_func, dist_type = dist_type)
#' 
#' @export

new_BayesMixture <- function(mcmc,
                             data,
                             K,
                             burnin,
                             dist = "NA",
                             pdf_func = NULL,
                             dist_type,
                             loglik = NULL) {
  ## input checks
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.string(dist_type),
              msg = "dist_type should be a string")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  assert_that(is.scalar(burnin), msg = "nb_iter should be an integer positive or zero")
  assert_that(burnin < nrow(mcmc),
              msg = "burnin parameter should be less than the number of mcmc draws")
  ##
  
  BayesMix = list(data = data,
                  dist_type = dist_type,
                  loglik = loglik)
  
  mcmc = as_draws_matrix(mcmc)

  # check that pars_names and mcmc match
  pars_names = unique(str_extract(str_to_lower(colnames(mcmc)), "[a-z]+"))

  if (dist == "poisson"){
    assert_that(sum(pars_names %in% c("eta", "lambda"))==2,
                msg = "variable names in mcmc output should be eta and lambda when dist = poisson")
  }
  
  if (dist == "shifted_poisson"){
    assert_that(sum(pars_names %in% c("eta", "kappa", "lambda"))==3,
                msg = "variable names in mcmc output should be eta and lambda when dist = shifted_poisson")
  }
  
  if (dist == "normal"){
    assert_that(sum(pars_names %in% c("eta", "mu", "sigma"))==3,
                msg = "variable names in mcmc output should be eta, mu and sigma when dist = normal")
  }
  
  if (dist == "skew_normal"){
    assert_that(sum(pars_names %in% c("eta", "xi", "omega", "alpha"))==4,
                msg = "variable names in mcmc output should be eta, xi, omega and alpha when dist = skew_normal")
  }
  
  ### arrange the mcmc matrix by variable type (mu1,mu2,...,muN,sigma...)
  # and select only the variables of interest
  mcmc_new = matrix(NA, nrow = nrow(mcmc), ncol = length(pars_names)*K)
  colnames(mcmc_new) = 1:ncol(mcmc_new)
  k_end = K
  k_start = 1
  for (par in pars_names) {
    
    #reorder
    cols_par = str_locate(colnames(mcmc), par)[,1]
    cols_par = which(!is.na(cols_par))
    mcmc_new[,k_start:k_end] = mcmc[, cols_par]
    
    #rename
    numb = gregexpr('[0-9]+', colnames(mcmc)[cols_par])
    numb = unlist(regmatches(colnames(mcmc)[cols_par],numb))
    
    colnames(mcmc_new)[k_start:k_end] = paste0(par, numb)
    
    k_start = k_end + 1
    k_end = k_start + K - 1
  }
  
  mcmc_all = mcmc_new
  mcmc = mcmc_all[(burnin+1):nrow(mcmc_all), ,drop = FALSE]
  
  
  if (dist %in% c("normal", "skew_normal",
                  "poisson", "shifted_poisson") & is.null(pdf_func)) {
    BayesMix$dist = dist
  } else {
    BayesMix$dist = "NA"
    BayesMix$pdf_func = pdf_func
    assert_that(!is.null(pdf_func),
                msg = "pdf_func is missing")
  }
  
  BayesMix$pars_names = pars_names
  BayesMix$mcmc = mcmc
  BayesMix$mcmc_all = mcmc_all
  
  class(BayesMix) <- "BayesMixture"
  
  return(BayesMix)
}