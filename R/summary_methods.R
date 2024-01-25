#' Summary method for \code{BayesMode} objects
#' 
#' @param object An object of class \code{BayesMode}.
#' @param ... Not used.
#' 
#' @export
summary.BayesMode <- function(object, ...) {
  modes = object$modes
  
  p1 = object$p1
  cat("The posterior probability of the data being multimodal is", 1-p1)
  
  tb_nb_modes = t(object$tb_nb_modes)
  colnames(tb_nb_modes) = c("Number of modes", "Posterior probabilty")
  tb_nb_modes = tb_nb_modes[order(tb_nb_modes[, 1]), ]

  cat("\n\n Number of estimated modes and their posterior probabilities:\n")
  tb_nb_modes
}


#' Summary method for \code{Mode} objects
#' 
#' @param object An object of class \code{Mode}.
#' @param ... Not used.
#' 
#' @export
summary.Mode <- function(object, ...) {
  Nb_m = length(object$mode_estimates)
  algo = object$algo
  d = object$dist
  K = object$K
  
  if (is.na(d)) {
    d = object$dist_type
  }
  
  if (Nb_m == 1) {
    m = "Mode"
  } else {
    m = "Modes"
  }
  
  cat("\n",m, "of a", d, "mixture with", K, "components.")
  cat("\n- Number of modes found:", Nb_m)
  cat("\n- Mode estimation technique:", object$algo, "algorithm")
}


#' Summary method for \code{Mixture} objects
#' 
#' @param object An object of class \code{Mode}.
#' @param ... Not used.
#' 
#' @export
summary.Mixture <- function(object, ...) {
  cat("\n Estimated mixture distribution.")
  cat("\n- Mixture type:", object$dist_type)
  cat("\n- Number of components:", object$K)
  cat("\n- Distribution family:", object$dist)
  cat("\n- Number of distribution variables:", object$nb_var)
  cat("\n- Names of variables:",
      object$pars_names[object$pars_names!="eta"])
}



#' Summary method for \code{BayesMixture} objects
#' The summary of MCMC draws is given by the function
#' \code{summarise_draws} from package \pkg{posterior}.
#' @param object An object of class \code{BayesMixture}.
#' @param ... Not used.
#' 
#' @importFrom posterior summarise_draws
#' 
#' @export
summary.BayesMixture <- function(object, ...) {
  d = object$dist
  K = object$K
  
  if (is.na(d)) {
    d = object$dist_type
  }
  
  cat("\n Mixture estimated with a Bayesian MCMC method.")
  cat("\n- Mixture type:", object$dist_type)
  cat("\n- Number of components:", object$K)
  cat("\n- Distribution family:", object$dist)
  cat("\n- Number of distribution variables:", object$nb_var)
  cat("\n- Names of variables:",
      object$pars_names[object$pars_names!="eta"])
  
  cat("\n\nSummary of MCMC output after burnin:\n")
  summarise_draws(object$mcmc)
}
