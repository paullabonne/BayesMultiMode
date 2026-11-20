#' Trace plots
#' 
#' This is wrapper around the [bayesplot::mcmc_trace()] function from package `bayesplot`.
#'
#' @param BayesMix An object of class `bayes_mixture`.
#' @param mcmc_vars Variables to plot. For telescoping the default is `c("K", "Kplus", "concentration")`; component columns can also be selected. For SFM the default is all columns.
#' @param with_burnin Plot all draws ?
#' @param ... Additional arguments passed to function [bayesplot::mcmc_trace()].
#' 
#' @importFrom bayesplot mcmc_trace
#' @importFrom assertthat assert_that
#' 
#' @return A trace plot.
#' 
#' @examples
#' # Example with galaxy data ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = galaxy
#'
#' # estimation
#' bayesmix = bayes_fit(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'
#' # trace plot
#' # bayes_trace(bayesmix)
#' 
#' @export
#'
bayes_trace <- function(BayesMix,
                        mcmc_vars = NULL,
                        with_burnin = FALSE,
                        ...) {
  
  assert_that(inherits(BayesMix, "bayes_mixture"), msg = "BayesMix should be an object of class bayes_mixture")
  assert_that(is.logical(with_burnin), msg = "with_burnin should be either TRUE or FALSE")
  
  if (with_burnin){
    mcmc = BayesMix$mcmc_all
  } else {
    mcmc = BayesMix$mcmc
  }
  if (!is.null(BayesMix$state)) {
    state <- if (with_burnin) BayesMix$state_all else BayesMix$state
    if (is.null(mcmc_vars)) mcmc_vars <- colnames(state)
    mcmc <- cbind(mcmc, state)[, mcmc_vars, drop = FALSE]
  }
  
  if (is.null(mcmc_vars)) {
    mcmc_vars = colnames(mcmc) 
  }
  
  message(cat("\nNote that label-switching might occur in the MCMC draws because BayesMultiMode does not carry out post-processing.",
              "\nWhile label-switching does not affect mode inference it can affect diagnostic checks."))
  mcmc_trace(mcmc, pars = mcmc_vars, ...)
}
