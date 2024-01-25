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
#' @param vars_to_keep (optional) Character vector containing the names
#' of the variables to keep in mcmc.

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
#' mu_mat = matrix(rep(mu, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = T)
#'
#' sigma = c(1,2)
#' sigma_mat = matrix(rep(sigma, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = T)
#' 
#' nu = c(5,5)
#' nu_mat = matrix(rep(nu, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = T)
#' 
#' eta = c(0.8,0.2)
#' eta_mat = matrix(rep(eta[1], 100) + rnorm(100, 0, 0.05),
#'             ncol = 1)
#' eta_mat = cbind(eta_mat,1-eta_mat)
#' 
#' xi_mat = matrix(0,100,2)
#' 
#' dist_type = "continuous"
#' 
#' data = c(sn::rst(eta[1]*1000, mu[1], sigma[1], nu = nu[1]),
#'         sn::rst(eta[2]*1000, mu[2], sigma[2], nu = nu[2]))
#' 
#' fit = cbind(eta_mat, mu_mat, sigma_mat, nu_mat, xi_mat)
#' colnames(fit) = c("eta1", "eta2", "mu1", "mu2",
#'                   "sigma1", "sigma2", "nu1", "nu2", "xi1", "xi2")
#' pdf_func = function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' BM = new_BayesMixture(fit, data, K = 2, burnin = 50,
#' pdf_func = pdf_func, dist_type = dist_type, loc = "xi")
#' # plot(BM)
#' @export

new_BayesMixture <- function(mcmc,
                             data,
                             K,
                             burnin,
                             dist = NA_character_,
                             pdf_func = NULL,
                             dist_type = NA_character_,
                             loglik = NULL,
                             vars_to_keep = NA_character_,
                             loc = NA_character_) {
  ## input checks
  assert_that(is.matrix(mcmc))
  assert_that(is.string(dist))
  assert_that(is.string(dist_type))
  if(!is.na(dist_type)) {
    assert_that(dist_type %in% c("continuous", "discrete"),
                msg = "dist_type should be either continuous or discrete")
  }
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(is.scalar(K) & K > 0, msg = "K should be a positive integer")
  assert_that(is.scalar(burnin) & burnin >= 0, msg = "burnin should be an integer positive or zero")
  assert_that(burnin < nrow(mcmc),
              msg = "burnin parameter should be less than the number of mcmc draws")
  assert_that(!(is.na(dist) & is.null(pdf_func)),
              msg = "one of dist or pdf_func must be specified")
  ## input checks
  assert_that(is.character(vars_to_keep))
  
  ##
  rownames(mcmc) = NULL
  mcmc_all = mcmc
  mcmc = mcmc_all[(burnin+1):nrow(mcmc_all), ,drop = FALSE]
  
  # extract parameter names
  col_names = str_extract(colnames(mcmc), "[a-z]+")
  pars_names = unique(col_names)
  
  # check that eta is included
  assert_that("eta" %in% pars_names,
              msg = "mcmc should include a parameter named eta representing mixture proportions.")
  
  # keep only variables specify
  if (sum(!is.na(vars_to_keep))>0) {
    pars_names = pars_names[pars_names %in% vars_to_keep]
    mcmc = mcmc[ , col_names %in% pars_names, drop = F]
  }
  
  # count number of components
  match = T
  K_from_names = rep(NA_real_, length(pars_names))
  for (i in 1:length(pars_names)) {
    K_from_names[i] = sum(col_names==pars_names[i])
  }
  
  assert_that(sum(K_from_names != K) == 0,
              msg = "There is a least one variable in mcmc that has not K components")

  list_func = test_and_export(mcmc[1,], pdf_func, dist, pars_names, dist_type, par_type = "mcmc", loc)

  BayesMix = list(mcmc = mcmc,
                  data = data,
                  mcmc_all = mcmc_all,
                  dist_type = list_func$dist_type,
                  loglik = loglik,
                  dist = dist,
                  pdf_func = list_func$pdf_func,
                  pars_names = pars_names,
                  loc = list_func$loc,
                  nb_var = length(pars_names) - 1, #minus the shares
                  K = ncol(mcmc)/length(pars_names))
  
  class(BayesMix) <- "BayesMixture"
  
  return(BayesMix)
}