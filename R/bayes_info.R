#' Compute the Bayesian Information Criterion
#'
#' Returns the BIC for each draw as well as its posterior mean and standard deviation.
#'
#' @param BayesMix An object of class `bayes_mixture`.
#'
#' @return Returns a list with:
#'  \item{BIC}{Vector of BIC values for each MCMC draw.}
#'  \item{mean_BIC}{Mean BIC}
#'  \item{sd_BIC}{Standard Deviation of the BIC}
#'
#' @examples
#' # Example with galaxy data ================================================
#' set.seed(123)
#'
#' # retrieve galaxy data
#' y <- galaxy
#'
#' # estimation
#' bayesmix <- bayes_fit(
#'   data = y,
#'   K = 5, # not many to run the example rapidly
#'   dist = "normal",
#'   nb_iter = 500, # not many to run the example rapidly
#'   burnin = 100
#' )
#'
#' # BIC
#' bayes_info(bayesmix)
#'
#' @export
#'
bayes_info <- function(BayesMix) {
  assert_that(inherits(BayesMix, "bayes_mixture"), msg = "BayesMix should be an object of class bayes_mixture")

  # inputs
  ll <- BayesMix$loglik
  n <- length(BayesMix$data)
  k <- BayesMix$nb_var

  # outputs
  BIC <- -2 * ll + k * log(n)
  AIC <- -2 * ll + 2 * k

  IC = cbind(BIC, AIC)
  colnames(IC) = c("BIC", "AIC")

  mean_IC <- apply(IC, 2, mean)
  sd_IC <- apply(IC, 2, sd)

  return(list(IC = IC, mean_IC = mean_IC, sd_IC = sd_IC))
}
