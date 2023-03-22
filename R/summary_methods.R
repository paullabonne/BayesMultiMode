#' Summary of Bayesian mode estimates.
#' 
#' @param object An object of class BayesMode.
#' @param ... Not used.
#' 
#' @export
summary.BayesMode <- function(object, ...) {
  stopifnot(inherits(object, "BayesMode"))
  
  modes = object$modes
  
  p1 = object$p1
  cat("The posterior probability of the data being multimodal is", 1-p1, ".")
  
  tb_nb_modes = t(object$tb_nb_modes)
  colnames(tb_nb_modes) = c("Number of modes", "Posterior probabilty")
  tb_nb_modes = tb_nb_modes[order(tb_nb_modes[, 1]), ]

  cat("\n\nThe number of estimated modes and their posterior probabilities is:\n")
  tb_nb_modes
}