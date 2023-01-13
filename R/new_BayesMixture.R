#' Creating a S3 object of class `BayesMixture`.
#' This function is helpful for users who want to explore modes in MCMC draws which have not been
#' derived using the function `bayes_estimation()`.
#' 
#' @param fit (numeric matrix) MCMC draws.
#' @param data (string) Distribution family of the mixture components supported by
#' the package (e.g. "normal", "student", "skew_normal", "shifted_poisson").
#' @param dist ...
#' @param pars_name (named character vector) providing the mapping between the distribution parameters names.
#' This input is used only if dist_name is invalid or NULL.
#' @param density (function) Pdf or pmf of the mixture components.
#' This input is used only if dist_name is invalid or NULL.
#' @param dist_type ...
#' 
#' @return An object of class `BayesMixture`.
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom stringr str_replace
#' 
#' @export

new_BayesMixture <- function(fit, data, dist = "NA", pars_name = "NA", density = NULL, dist_type = "NA") {
  stopifnot(is.character(pars_name))
  
  BayesMix = list(data = data,
                  fit = fit,
                  dist_type = dist_type)
  
  mcmc <- as_draws_matrix(fit)
  
  if (dist == "shifted_poisson") {
    mcmc = post_sfm_mcmc(mcmc)
  }
  
  if (length(pars_name) > 1) {
    # change variable names
    for (i in 1:length(pars_name)) {
      colnames(mcmc) = str_replace(colnames(mcmc), names(pars_name)[i], pars_name[i])
    }
  }

  if (dist %in% c("normal", "student", "skew_normal", "shifted_poisson")) {
    BayesMix$dist = dist
  } else {
    BayesMix$dist = "NA"
    BayesMix$density = density
  }
  
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