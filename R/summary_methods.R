#' Summary method for `bayes_mode` objects
#'
#' @param object An object of class `bayes_mode`.
#' @param ... Not used.
#'
#' @export
summary.bayes_mode <- function(object, ...) {
  modes <- object$modes

  if (!is.null(object$conditional_nb_modes)) {
    cat("\n", "These results are conditional on the number of modes being", object$conditional_nb_modes, "\n")
  } else {
    p1 <- object$p1
    cat("The posterior probability of multimodality is", 1 - p1, "\n")

    cat(
      "\n The most likely number of modes is",
      t(object$p_nb_modes)[t(object$p_nb_modes)[, 2] == max(t(object$p_nb_modes)[, 2]), 1],
      "\n"
    )

    cat("\nInference results on the number of modes:")
    cat("\n  p_nb_modes")
    head_print(t(object$p_nb_modes))
  }

  cat("\nInference results on mode locations:")
  cat("\n  p_loc")
  head_print(t(object$p_mode_loc))
}


#' Summary method for `mix_mode` objects
#'
#' @param object An object of class `mix_mode`.
#' @param ... Not used.
#'
#' @export
summary.mix_mode <- function(object, ...) {
  Nb_m <- length(object$mode_estimates)
  algo <- object$algo
  d <- object$dist
  K <- object$K

  if (is.na(d)) {
    d <- object$dist_type
  }

  if (Nb_m == 1) {
    m <- "Mode"
  } else {
    m <- "Modes"
  }

  cat(m, "of a", d, "mixture with", K, "components.")
  cat("\n- Number of modes found:", Nb_m)
  cat("\n- Mode estimation technique:", object$algo, "algorithm")
  cat("\n- Estimates of mode locations:")
  cat("\n  mode_estimates")
  head_print(round(object$mode_estimates), 3)
}


#' Summary method for `mixture` objects
#'
#' @param object An object of class `mixture`.
#' @param ... Not used.
#'
#' @export
summary.mixture <- function(object, ...) {
  cat("Estimated mixture distribution.")
  cat("\n- Mixture type:", object$dist_type)
  cat("\n- Number of components:", object$K)
  cat("\n- Distribution family:", object$dist)
  cat("\n- Number of distribution variables:", object$nb_var)
  cat(
    "\n- Names of variables:",
    object$pars_names[object$pars_names != "eta"]
  )
  cat("\n- Parameter estimates:")
  cat("\n  pars")
  head_print(object$pars)
}



#' Summary method for `bayes_mixture` objects
#' The summary of MCMC draws is given by the function
#' `summarise_draws` from package \pkg{posterior}.
#' @param object An object of class `bayes_mixture`.
#' @param ... Not used.
#'
#' @importFrom posterior summarise_draws
#'
#' @export
summary.bayes_mixture <- function(object, ...) {
  d <- object$dist
  K <- object$K

  if (is.na(d)) {
    d <- object$dist_type
  }

  cat("Mixture estimated with a Bayesian MCMC method.")
  cat("\n- Mixture type:", object$dist_type)
  cat("\n- Number of components:", object$K)
  cat("\n- Distribution family:", object$dist)
  cat("\n- Number of distribution variables:", object$nb_var)
  cat(
    "\n- Names of variables:",
    object$pars_names[object$pars_names != "eta"]
  )

  cat("\n\nSummary of MCMC output after burnin:\n")
  print(summarise_draws(object$mcmc))
  cat(paste0("this table can be reproduced with: summarise_draws(", deparse(substitute(object)), "$mcmc)"))
  message(cat(
    "\n\nNote that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing.",
    "\nWhile label-switching does not affect mode inference it can affect diagnostic checks.\n"
  ))
}
