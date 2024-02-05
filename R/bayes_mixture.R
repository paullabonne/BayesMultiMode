#' Creating a S3 object of class `bayes_mixture`
#' 
#' Creates an object of class `bayes_mixture` which can subsequently be used as argument in [bayes_mode()].
#' This function is useful for users who want to use the mode inference capabilities of `BayesMultiMode` with mixture
#' estimated using external software.
#' 
#' @param mcmc A matrix of MCMC draws with one column per variable, e.g. eta1, eta2, ..., mu1, mu2, etc...
#' @param data Vector of observation used for estimating the model.
#' @param burnin Number of draws to discard as burnin.
#' @param dist Distribution family of the mixture components supported by
#' the package (i.e. `"normal"`, `"student"`, `"skew_normal"` or `"shifted_poisson"`).
#' If left unspecified, `pdf_func` is required.
#' @param pdf_func (function) Pdf or pmf of the mixture components;
#' this input is used only if `dist` is left unspecified.
#' pdf_func should have two arguments : (i) the observation where the pdf is evaluated;
#' (ii) a named vector representing the function parameters. For instance a normal pdf would take the form:
#' `pdf_func <- function(x, pars) dnorm(x, pars['mu'], pars['sigma'])`.
#' The names of `pars` should correspond to variables in `mcmc`, e.g. `"mu1"`, `"mu2"` etc... 
#' @param dist_type Either `"continuous"` or `"discrete"`.
#' @param loglik Vector showing the log likelihood at each MCMC draw.
#' @param vars_to_keep (optional) Character vector containing the names
#' of the variables to keep in `mcmc`.
#' @param vars_to_rename (optional) Use for renaming variables/parameters in `mcmc`.
#' A named character vector where the names are the new variable names
#' and the elements the variables in `mcmc`, e.g. c("new_name" = "old_name").
#' @param loc (for continuous mixtures other than Normal mixtures) String indicating the location parameter
#' of the distribution; the latter is used to initialise the MEM algorithm.
#' 
#' @return A list of class `bayes_mixture` containing:
#'  \item{data}{Same as argument.}
#'  \item{mcmc}{Matrix of MCMC draws where the rows corresponding to burnin have been discarded;}
#'  \item{mcmc_all}{Matrix of MCMC draws.}
#'  \item{loglik}{Log likelihood at each MCMC draw.}
#'  \item{K}{Number of components.}
#'  \item{dist}{Same as argument.}
#'  \item{pdf_func}{The pdf/pmf of the mixture components.}
#'  \item{dist_type}{Type of the distribution, i.e. continuous or discrete.}
#'  \item{pars_names}{Names of the mixture components' parameters.}
#'  \item{loc}{Name of the location parameter of the mixture components.}
#'  \item{nb_var}{Number of parameters in the mixture distribution.}
#' 
#' @importFrom posterior subset_draws
#' @importFrom stringr str_extract
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_locate
#' 
#' @examples
#' 
#' # Example with a Student t ================================================
#' 
#' # Constructing synthetic mcmc output
#' mu = c(0.5,6)
#' mu_mat = matrix(rep(mu, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = TRUE)
#'
#' omega = c(1,2)
#' sigma_mat = matrix(rep(omega, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = TRUE)
#' 
#' nu = c(5,5)
#' nu_mat = matrix(rep(nu, 100) + rnorm(200, 0, 0.1),
#'             ncol = 2, byrow = TRUE)
#' 
#' eta = c(0.8,0.2)
#' eta_mat = matrix(rep(eta[1], 100) + rnorm(100, 0, 0.05),
#'             ncol = 1)
#' eta_mat = cbind(eta_mat,1-eta_mat)
#' 
#' xi_mat = matrix(0,100,2)
#' 
#' fit = cbind(eta_mat, mu_mat, sigma_mat, nu_mat, xi_mat)
#' colnames(fit) = c("eta1", "eta2", "mu1", "mu2",
#'                   "omega1", "omega2", "nu1", "nu2", "xi1", "xi2")
#'                   
#' # sampling observations
#' data = c(sn::rst(eta[1]*1000, mu[1], omega[1], nu = nu[1]),
#'         sn::rst(eta[2]*1000, mu[2], omega[2], nu = nu[2]))
#'         
#' pdf_func = function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' dist_type = "continuous"
#' 
#' BM = bayes_mixture(fit, data, burnin = 50,
#' pdf_func = pdf_func, dist_type = dist_type,
#' vars_to_rename = c("sigma" = "omega"), loc = "xi")
#' # plot(BM)
#' @export

bayes_mixture <- function(mcmc,
                          data,
                          burnin,
                          dist = NA_character_,
                          pdf_func = NULL,
                          dist_type = NA_character_,
                          loglik = NULL,
                          vars_to_keep = NA_character_,
                          vars_to_rename = NA_character_,
                          loc = NA_character_) {
  ## input checks
  assert_that(is.matrix(mcmc))
  assert_that(is.string(dist))
  assert_that(is.string(dist_type))
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(is.scalar(burnin) & burnin >= 0, msg = "burnin should be an integer positive or zero")
  assert_that(burnin < nrow(mcmc),
              msg = "burnin parameter should be less than the number of mcmc draws")
  ## input checks
  assert_that(is.character(vars_to_keep))
  assert_that(is.character(vars_to_rename))
  ##
  rownames(mcmc) = NULL
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
  
  if(any(!is.na(vars_to_rename))) {
    assert_that(!is.null(names(vars_to_rename)),
                msg = "vars_to_rename should be named character vector")
    assert_that(all(vars_to_rename %in% pars_names),
                msg = "old variable names in vars_to_rename should all be in the retained mcmc variables")
    
    new_names = colnames(mcmc)
    for (i in 1:length(vars_to_rename)) {
      new_names = str_replace_all(new_names,
                                  vars_to_rename[[i]],
                                  names(vars_to_rename)[i])
    }
    
    colnames(mcmc) = new_names
    pars_names = unique(str_extract(new_names, "[a-z]+"))
  }
  
  list_func = test_and_export(mcmc[1,,drop =T], pdf_func, dist, pars_names, dist_type, loc)
  
  BayesMix = list(data = data,
                  mcmc = mcmc,
                  mcmc_all = mcmc_all,
                  loglik = loglik,
                  K = list_func$K,
                  dist = dist,
                  dist_type = list_func$dist_type,
                  pdf_func = list_func$pdf_func,
                  pars_names = pars_names,
                  loc = list_func$loc,
                  nb_var = length(pars_names) - 1) #minus the shares
  
  class(BayesMix) <- "bayes_mixture"
  
  return(BayesMix)
}