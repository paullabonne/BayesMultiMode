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
#' @param pars_names Names of the mixture parameters; first element should 
#' correspond to the mixture proportions.
#' @param pdf_func Pdf or pmf of the mixture components;
#' this input is used only if dist_name is invalid or NULL.
#' @param dist_type Either "continuous" or "discrete".
#' 
#' @returns
#' A list of class \code{BayesMixture} containing:
#' \itemize{
#'  \item{data}{ - Same as argument.}
#'  \item{dist_type}{ - Same as argument.}
#'  \item{pars_names}{ - Same as argument.}
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
#' BM = new_BayesMixture(fit, data, K = 2, burnin = 1,
#' pars_names = pars_names, pdf_func = pdf_func, dist_type = dist_type)
#' 
#' @export

new_BayesMixture <- function(mcmc,
                             data,
                             K,
                             burnin,
                             dist = "NA",
                             pars_names,
                             pdf_func = NULL,
                             dist_type) {
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
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##
  
  BayesMix = list(data = data,
                  dist_type = dist_type)
  
  mcmc = as_draws_matrix(mcmc)

  # check that pars_names and mcmc match
  names_mcmc = str_to_lower(colnames(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")

  if (dist == "NA") {
    
    assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
                msg = "the name of the parameters provided by pars_names and those of the mcmc object do not match")
    
  } else {
    
    if (dist == "poisson" & sum(c("eta", "lambda") %in% names_mcmc) < 2){
      assert_that(sum(pars_names %in% c("eta", "lambda"))==2,
                  msg = "the name of the parameters provided by pars_names should be eta and lambda")
      assert_that(sum(pars_names %in% names_mcmc)==2,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    if (dist == "normal" & sum(c("eta", "mu", "sigma") %in% names_mcmc)<3){
      assert_that(sum(pars_names %in% c("eta", "mu", "sigma"))==3,
                  msg = "the name of the parameters provided by pars_names should be eta, mu and sigma")
      assert_that(sum(pars_names %in% names_mcmc)==3,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    if (dist == "skew_normal" & sum(c("eta", "xi", "omega", "alpha") %in% names_mcmc)<4){
      assert_that(sum(pars_names %in% c("eta", "mu", "sigma", "xi"))==4,
                  msg = "the name of the parameters provided by pars_names should be eta, mu, sigma and xi")
      assert_that(sum(pars_names %in% names_mcmc)==4,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
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