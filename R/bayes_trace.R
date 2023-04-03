#' Trace plots
#' 
#' This is wrapper around the `mcmc_trace()` function from package bayesplot.
#'
#' @param BayesMix An object of class BayesMixture
#' @param mcmc_vars Variables to plot; default is all the variable in the MCMC output
#' @param with_burnin Plot all draws ?
#' @param ... Additional arguments passed to function `mcmc_trace()` from the package bayesplot.
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
#' bayesmix = bayes_estimation(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'
#' # trace plot
#' bayes_trace(bayesmix)
#' 
#' @export
#'
bayes_trace <- function(BayesMix,
                        mcmc_vars = NULL,
                        with_burnin = FALSE,
                        ...) {
  
  assert_that(inherits(BayesMix, "BayesMixture"), msg = "BayesMix should be an object of class BayesMixture")
  assert_that(is.logical(with_burnin), msg = "with_burnin should be either TRUE or FALSE")
  
  if (with_burnin){
    mcmc = BayesMix$mcmc_all
  } else {
    mcmc = BayesMix$mcmc
  }
  
  if (is.null(mcmc_vars)) {
    mcmc_vars = colnames(mcmc) 
  }
  
  mcmc_trace(mcmc, pars = mcmc_vars, ...)
}
