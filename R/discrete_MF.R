#' Mode-finding algorithm for mixture of discrete distributions
#' 
#' Function to estimate modes in mixtures of discrete distributions following; see Cross et al. (2023).
#' 
#' @param mixture An object of class Mixture.
#' @param type Type of modes, either unique or all (the latter includes flat modes); default is "all".
#' 
#' @return An object of class Mode.
#' 
#' @references
#' \insertRef{cross_2023}{BayesMultiMode}
#' 
#' @details
#' 
#' This algorithm returns the local maxima of the mixture
#' \deqn{p(x) = \sum_{k=1}^{K}\pi_k p_k(x).}
#' 
#'  By definition, modes must satisfy either: 
#'  \deqn{p_k(y_{m}-1) < p_k(y_{m}) > p_k(y_{m}+1)};
#'  \deqn{p_k(y_{m}-1) < p_k(y_{m}) = p_k(y_{m}+1) = \ldots = p_k(y_{m}+l-1) > p_k(y_{m}+l).}
#'  
#'  The algorithm evaluate each location point with these two conditions.
#' 
#' @references
#' \insertRef{schaap_genome-wide_2013}{BayesMultiMode}
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
#' 
#' @examples
#' # Example with the poisson distribution ====================================
#' lambda = c(0.1,10)
#' p = c(0.5,0.5)
#' params = c(eta = p, lambda = lambda)
#' dist = "poisson"
#' 
#' data = c(rpois(p[1]*1e3, lambda[1]),
#'          rpois(p[2]*1e3, lambda[2]))
#' 
#' mix = new_Mixture(params, data = data, dist = dist)
#' 
#' modes = discrete_MF(mix)
#'
#' # summary(modes)
#' # plot(modes, from = 0, to = 20)
#' 
#' # Example with an arbitrary distribution ===================================
#' mu = c(20,5)
#' size = c(20,0.5)
#' p = c(0.5,0.5)
#' params = c(eta = p, mu = mu, size = size)
#' 
#' data = c(rnbinom(p[1]*1e3, mu = mu[1], size = size[1]),
#'          rnbinom(p[2]*1e3, mu = mu[2], size = size[2]))
#' 
#' pmf_func <- function(x, pars) {
#'   dnbinom(x, mu = pars["mu"], size = pars["size"])
#' }
#' 
#' mix = new_Mixture(params, data = data, pdf_func = pmf_func, dist_type = "discrete")
#' modes = discrete_MF(mix)
#' 
#' # summary(modes)
#' # plot(modes, from = 0, to = 50)
#' 
#' @export

discrete_MF <- function(mixture, type = "all"){
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  dist = mixture$dist
  pdf_func = mixture$pdf_func
  data = mixture$data
  
  ## input checks
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(!any(is.na(data)) & !any(is.infinite(data)),
              msg = "data should not include missing or infinite values")
  assert_that(type %in% c("unique", "all"),
              msg = "type must be either 'unique' or 'all' ")
  ##
  
  ##
  x = min(data):max(data)
  ##
  
  pars_mat <- vec_to_mat(pars, pars_names)
  ##
  
  Khat = nrow(pars_mat)
  
  ### Getting denisty
  py = pdf_func_mix(x, pars_mat, pdf_func)

  # change in the pdf
  d_py = diff(py)
  
  # where does the pdf decrease ?
  x_decrease = x[d_py<0]

  # Only keep the points where the pdf starts to decrease; these are modes
  d2_py = c(0, x_decrease[-1] - x_decrease[-length(x_decrease)])
  x_decrease = x_decrease[which(d2_py!=1)]
  
  # get pdf at these modes 
  pdf_modes = py[x %in% x_decrease]
  
  # get all the points at these peaks (there might be flat modes) 
  loc_modes = x[which(py %in% pdf_modes)]
  
  if (length(loc_modes) != length(x_decrease)) {
    warning("Some modes are flat.")
  }
  
  output = rep(NA_real_, length(x))
  
  if (type == "unique") {
    output[1:length(loc_modes)] = x_decrease
  }
  
  if (type == "all") {
    output[1:length(loc_modes)] = loc_modes
  }
  
  mode = list()
  mode$mode_estimates = output[!is.na(output)]
  mode$dist = dist
  mode$pars = pars
  mode$pdf_func = pdf_func
  mode$data = data
  mode$dist_type = "discrete"
  mode$py = py
  mode$algo = "discrete"
  mode$K = mixture$K
  mode$nb_var = mixture$nb_var
  
  class(mode) = "Mode"

  return(mode)
}