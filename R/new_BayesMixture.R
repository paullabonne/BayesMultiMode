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
#' fit = c(eta = p, mu = mu, sigma = sigma, nu = nu, xi = c(0,0))
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
                             dist = NA_character_,
                             pdf_func = NULL,
                             dist_type = NA_character_,
                             loglik = NULL,
                             vars_to_keep = NA_character_) {
  ## input checks
  assert_that(is.matrix(mcmc),
              msg = "new_BayesMixture failed; mcmc should be a matrix")
  assert_that(is.string(dist),
              msg = "new_BayesMixture failed; dist should be a string")
  assert_that(is.string(dist_type),
              msg = "new_BayesMixture failed; dist_type should be a string")
  if(!is.na(dist_type)) {
    assert_that(dist_type %in% c("continuous", "discrete"),
                msg = "dist_type should be either continuous or discrete")
  }
  assert_that(is.vector(data) & length(data) > 0,
              msg = "new_BayesMixture failed; data should be a vector of length > 0")
  assert_that(is.scalar(K) & K > 0, msg = "new_BayesMixture failed; K should be a positive integer")
  assert_that(is.scalar(burnin), msg = "new_BayesMixture failed; nb_iter should be an integer positive or zero")
  assert_that(burnin < nrow(mcmc),
              msg = "new_BayesMixture failed; burnin parameter should be less than the number of mcmc draws")
  assert_that(!(is.na(dist) & is.null(pdf_func)),
              msg = "you have to specify either dist or pdf_func")
  ## input checks
  assert_that(is.character(vars_to_keep))
  
  ##
  mcmc_all = mcmc
  mcmc = mcmc_all[(burnin+1):nrow(mcmc_all), ,drop = FALSE]
  
  # extract parameter names
  col_names = str_extract(colnames(mcmc), "[a-z]+")
  pars_names = unique(col_names)
  
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
              msg = "new_BayesMixture failed; There is a least one variable in mcmc that has not K components")
  
  if (!is.na(dist)) {
    if (dist == "poisson"){
      assert_that(sum(pars_names %in% c("eta", "lambda"))==2,
                  msg = "new_BayesMixture failed; variable names in mcmc output should be eta and lambda when dist = poisson")
    }
    
    if (dist == "shifted_poisson"){
      assert_that(sum(pars_names %in% c("eta", "kappa", "lambda"))==3,
                  msg = "new_BayesMixture failed; variable names in mcmc output should be eta and lambda when dist = shifted_poisson")
    }
    
    if (dist == "normal"){
      assert_that(sum(pars_names %in% c("eta", "mu", "sigma"))==3,
                  msg = "new_BayesMixture failed; variable names in mcmc output should be eta, mu and sigma when dist = normal")
    }
    
    if (dist == "skew_normal"){
      assert_that(sum(pars_names %in% c("eta", "xi", "omega", "alpha"))==4,
                  msg = "new_BayesMixture failed; variable names in mcmc output should be eta, xi, omega and alpha when dist = skew_normal")
    }
    
    if (dist %in% c("normal", "skew_normal")) {
      dist_type = "continuous"
    } else if (dist %in% c("poisson", "shifted_poisson")) {
      dist_type = "discrete"
    } else {
      stop("Unsupported distribution; dist should be either normal, skew_normal, poisson or shifted_poisson")
    } 
  }
  
  # that pdf_func can be computed when provided
  if (!is.null(pdf_func)) {
    assert_that(!is.na(pdf_func(1, vec_to_mat(mcmc[1, ], pars_names)[1,-1])),
                msg = "new_BayesMixture failed; running pdf_func with pars provided returns NA") 
    assert_that(!is.na(dist_type),
                msg = "new_BayesMixture failed; dist_type must be provided when argument pdf_func is used") 
  }
  
  BayesMix = list(mcmc = mcmc,
                  data = data,
                  mcmc_all = mcmc_all,
                  dist_type = dist_type,
                  loglik = loglik,
                  dist = dist,
                  pdf_func = pdf_func,
                  pars_names = pars_names)
  
  class(BayesMix) <- "BayesMixture"
  
  return(BayesMix)
}