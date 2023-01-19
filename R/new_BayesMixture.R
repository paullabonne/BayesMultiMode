#' Creating a S3 object of class `BayesMixture`.
#' This function is helpful for users who want to explore modes in MCMC draws which have not been
#' derived using the function `bayes_estimation()`.
#' 
#' @param fit A matrix of MCMC draws.
#' @param data A vector containing the data used for estimating the model and generating the MCMC draws.
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
#' 
#' @export

new_BayesMixture <- function(fit, data, dist = "NA", pars_names, pdf_func = NULL, dist_type) {
  ## input checks
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.string(dist_type),
              msg = "dist_type should be a string")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  ##
  
  BayesMix = list(data = data,
                  fit = fit,
                  dist_type = dist_type)
  
  mcmc <- as_draws_matrix(fit)
  
  if (dist == "shifted_poisson") {
    mcmc_par <- attributes(mcmc)
    #Burn in
    mcmc = mcmc[(mcmc_par$warmup+1):nrow(mcmc), ]
  }
  
  # check that pars_names and mcmc match
  names_mcmc = str_to_lower(colnames(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")

  # adapt cols names
  col_names = names_mcmc
  j = 1
  name_previous = "NA"
  
  for (i in 1:length(col_names)) {
    if (col_names[i] == name_previous) {
      j = j+1
    } else {
      j = 1
    }
    name_previous = col_names[i]
    col_names[i] = paste0(col_names[i], "[", j, "]")
  }
  colnames(mcmc) <- col_names
  
  names_mcmc = unique(names_mcmc)
  
  if (dist == "NA") {
    assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
                msg = "the name of the parameters provided by pars_names and those of the mcmc object do not match") 
  } else {
    change = FALSE
    
    if (dist == "normal" & sum(c("theta", "mu", "sigma") %in% names_mcmc)<3){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma"))==3,
                  msg = "the name of the parameters provided by pars_names should be theta, mu and sigma")
      assert_that(sum(names_mcmc %in% pars_names)==3,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    if (dist == "student" & sum(c("theta", "mu", "sigma", "nu") %in% names_mcmc)<4){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma and nu")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma", "nu"))==4,
                  msg = "the name of the parameters provided by pars_names should be theta, mu, sigma and nu")
      assert_that(sum(names_mcmc %in% pars_names)==4,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    if (dist == "skew_normal" & sum(c("theta", "mu", "sigma", "xi") %in% names_mcmc)<4){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma and xi")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma", "xi"))==4,
                  msg = "the name of the parameters provided by pars_names should be theta, mu, sigma and xi")
      assert_that(sum(names_mcmc %in% pars_names)==4,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    if (dist == "skew_t" & sum(c("theta", "mu", "sigma", "xi", "nu") %in% names_mcmc)<5){
      change = TRUE
      
      assert_that(!is.null(names(pars_names)),
                  msg = "pars_names should be a named vector with names : theta, mu, sigma, xi and nu")
      assert_that(sum(names(pars_names) %in% c("theta", "mu", "sigma", "xi", "nu"))==5,
                  msg = "the name of the parameters provided by pars_names should be theta, mu, sigma, xi and nu")
      assert_that(sum(names_mcmc %in% pars_names)==5,
                  msg = "the name of the parameters provided by pars_names do match with the mcmc parameters") 
    }
    
    # keep only relevant variables
    mcmc = subset_draws(mcmc, variable = pars_names)
    
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
  
  if (dist %in% c("normal", "student", "skew_normal", "skew_t", "shifted_poisson")) {
    BayesMix$dist = dist
  } else {
    BayesMix$dist = "NA"
    BayesMix$pdf_func = pdf_func
  }
  
  BayesMix$pars_names = pars_names
  BayesMix$mcmc = mcmc
  
  class(BayesMix) <- "BayesMixture"
  
  return(BayesMix)
}