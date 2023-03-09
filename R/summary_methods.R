#' Summary modes
#' @param x ...
#' @param ... ...
#' 
#' @export
summary.BayesMode <- function(x, ...) {
  stopifnot(inherits(x, "BayesMode"))
  
  modes = x$modes
  
  p1 = x$p1
  cat("The posterior probability of the data being multimodal is", 1-p1, ".")
  
  tb_nb_modes = t(x$tb_nb_modes)
  colnames(tb_nb_modes) = c("Number of modes", "Posterior probabilty")
  tb_nb_modes = tb_nb_modes[order(tb_nb_modes[, 1]), ]

  cat("\n\nThe number of estimated modes and their posterior probabilities is:\n")
  show(tb_nb_modes)
}