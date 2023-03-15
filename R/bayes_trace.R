#' Trace plots. This is wrapper around the mcmc_trace function from package bayesplot.
#'
#' @param BayesMix An object of class BayesMixture.
#' @param mcmc_vars Variables to plot. Default is all the variable in the MCMC output.
#' @param ... arguments passed to mcmc_trace.
#' 
#' @importFrom bayesplot mcmc_trace
#' @importFrom assertthat assert_that
#' 
#' @return A trace plot.
#' 
#' @export
#'
#'
bayes_trace <- function(BayesMix,
                        mcmc_vars = NULL,
                        ...) {
  
  assert_that(inherits(BayesMix, "BayesMixture"), msg = "BayesMix should be an object of class BayesMixture")
  
  if (is.null(mcmc_vars)) {
    mcmc_vars = colnames(BayesMix$mcmc) 
  }
  
  mcmc_trace(BayesMix$mcmc, pars = mcmc_vars)
}
