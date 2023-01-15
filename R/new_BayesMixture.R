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
    mcmc = post_sfm_mcmc(mcmc)
  }
  
  
  # check that pars_names and mcmc match
  names_mcmc = str_to_lower(colnames(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
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
  
  if (dist %in% c("normal", "student", "skew_normal", "shifted_poisson")) {
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


#' @keywords internal
post_sfm_mcmc <- function(mcmc){
  
  mcmc_post = mcmc
  
  #Number of draws
  M = nrow(mcmc)
  
  mcmc_par <- attributes(mcmc)
  #Burn in
  mcmc_post = mcmc_post[(mcmc_par$warmup+1):M,]
  
  # Discard empty components
  # when a component is empty in a given draw it has a NA. 
  # Thus we discard components we are always NAs (empty in all draws).
  temp = mcmc_post
  temp[!is.na(temp)] = 1
  temp[is.na(temp)] = 0
  sum_temp = apply(temp,2,sum)
  
  mcmc_post = mcmc_post[,sum_temp>0]
  
  # number of non-tempty components
  # if(sfm_mcmc$mixt=="shifted_poisson"){
  #   K_non_emp = length(which(sum_temp>0))/3
  # }
  
  return(mcmc_post)
}