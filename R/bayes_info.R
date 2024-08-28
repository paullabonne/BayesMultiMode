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

  mcmc <- BayesMix$mcmc

  ll_full <- 0 # Initialize log_predictive_density to zero
  nb_draw <- nrow(mcmc) # Number of MCMC draws
  data = BayesMix$data
  n <- length(data) # Number of data points
  K <- BayesMix$K # Number of components in the mixture
  pdf_func = BayesMix$pdf_func
  pars_names = BayesMix$pars_names

  ll_full = rep(NA, nb_draw)

  # Loop over each MCMC sample
  for (i in 1:nb_draw) {
    # Get the parameters
    pars_i <- mcmc[i, ] # parameters

    pars_mat = vec_to_mat(pars_i, pars_names)
    pars_mat = na.omit(pars_mat)

    density_i = pdf_func_mix(data, pars_mat, pdf_func)

    density_i[density_i == 0] = 1e-100

    ######
    ll_full[i] = sum(log(density_i))
    ######
  }

  # Average over all MCMC draws
  mean_mll <- mean(ll_full)
  sd_mll <- sd(ll_full)
  return(list(mean_mll = mean_mll, sd_mll = sd_mll))
}
