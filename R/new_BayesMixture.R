#' Creating a S3 object of class `BayesMixture`.
#' This function is helpful for users who want to explore modes in MCMC draws which have not been
#' derived using the function `bayes_estimation()`.
#' 
#' @param mcmc (numeric matrix) MCMC draws.
#' @param dist_name (string) Distribution family of the mixture components supported by
#' the package (e.g. "normal", "student", "skew_normal", "shifted_poisson").
#' @param pars_name (named character vector) providing the mapping between the distribution parameters names.
#' This input is used only if dist_name is invalid or NULL.
#' @param density (function) Pdf or pmf of the mixture components.
#' This input is used only if dist_name is invalid or NULL.
#' 
#' @return An object of class `BayesMixture`.
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom stringr str_replace

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
