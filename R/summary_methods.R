#' Summary method for \code{BayesMode} objects
#' 
#' @param object An object of class \code{BayesMode}.
#' @param ... Not used.
#' 
#' @export
summary.BayesMode <- function(object, ...) {
  stopifnot(inherits(object, "BayesMode"))
  
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
  stopifnot(inherits(object, "Mode"))
  
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
  cat("\n", Nb_m, m, "of a", d,"distribution with", K, "components.")
  cat("\n Mode estimation using the", algo,"algorithm.")
}
