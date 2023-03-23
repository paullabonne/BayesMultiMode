#' Creating a S3 object of class `BayesMixture`.
#' 
#' This function is for users who want to explore modes in MCMC draws which have not been
#' derived using the function bayes_estimation().
#' 
#' @param mcmc A matrix of MCMC draws.
#' @param data A vector containing the data used for estimating the model and generating the MCMC draws.
#' @param K Number of mixture components.
#' @param burnin Number of draws to discard as burnin.
#' @param dist Distribution family of the mixture components supported by
#' the package (e.g. "normal", "student", "skew_normal", "shifted_poisson").
#' @param pars_names Mapping between the distribution parameters names.
#' This input is used only if dist_name is invalid or NULL.
#' @param pdf_func Pdf or pmf of the mixture components.
#' This input is used only if dist_name is invalid or NULL.
#' @param dist_type Either "continous" or "discrete"
#' 
#' @return An object of class `BayesMixture`.
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom posterior subset_draws
#' @importFrom stringr str_extract
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_replace
#' @importFrom stringr str_locate
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
    
    change = FALSE
    
    if (dist == "poisson" & sum(c("theta", "lambda") %in% names_mcmc) < 2){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta and lambda")
      assert_that(sum(names(pars_names) %in% c("theta", "lambda"))==2,
                  msg = "the name of the parameters provided by pars_names should be theta and lambda")
      assert_that(sum(pars_names %in% names_mcmc)==2,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    if (dist == "normal" & sum(c("theta", "mu", "sigma") %in% names_mcmc)<3){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma"))==3,
                  msg = "the name of the parameters provided by pars_names should be theta, mu and sigma")
      assert_that(sum(pars_names %in% names_mcmc)==3,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    if (dist == "skew_normal" & sum(c("theta", "xi", "omega", "alpha") %in% names_mcmc)<4){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma and xi")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma", "xi"))==4,
                  msg = "the name of the parameters provided by pars_names should be theta, mu, sigma and xi")
      assert_that(sum(pars_names %in% names_mcmc)==4,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    # change the parameter names in mcmc using pars_names
    if (change) {
      assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
                  msg = "the name of the parameters provided by pars_names and those of the mcmc object do not match") 
      
      # change variable names
      for (i in 1:length(pars_names)) {
        colnames(mcmc) = str_replace(colnames(mcmc), pars_names[i], names(pars_names)[i])
      }
    } 
  }
  
  ### arrange the mcmc matrix by variable type (mu1,mu2,...,muN,sigma...)
  # and select only the variables of interest
  mcmc_new = matrix(NA, nrow = nrow(mcmc), ncol = length(pars_names)*K)
  colnames(mcmc_new) = 1:ncol(mcmc_new)
  k_end = K
  k_start = 1
  for (par in names(pars_names)) {

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
  mcmc = mcmc_all[(burnin+1):nrow(mcmc_all), ]
  
  
  if (dist %in% c("normal", "skew_normal",
                  "poisson", "shifted_poisson") & is.null(pdf_func)) {
    BayesMix$dist = dist
  } else {
    BayesMix$dist = "NA"
    BayesMix$pdf_func = pdf_func
  }
  
  BayesMix$pars_names = pars_names
  BayesMix$mcmc = mcmc
  BayesMix$mcmc_all = mcmc_all
  
  class(BayesMix) <- "BayesMixture"
  
  return(BayesMix)
}