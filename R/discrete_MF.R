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
#' modes = discrete_MF(mix)
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
#' mix = new_Mixture(params, data = data, pdf_func = pmf_func)
#' modes = discrete_MF(mix)
#' 
#' @export

discrete_MF <- function(mixture, type = "all"){
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  dist = mixture$dist
  pmf_func = mixture$pdf_func
  data = mixture$data
  
  ## input checks
  assert_that(is.vector(pars),
              msg = "pars should be a vector")
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(!any(is.na(data)) & !any(is.infinite(data)),
              msg = "y should not include missing or infinite values")
  assert_that(type %in% c("unique", "all"),
              msg = "type should be either 'unique' or 'all' ")
  ##
  
  ##
  x = min(data):max(data)
  ##
  
  pars_mat <- vec_to_mat(pars, pars_names)
  ##
  
  Khat = nrow(pars_mat)
  
  ### Getting individual component densities
  
  pdf_k = matrix(0, nrow=length(x), ncol=Khat) 
  for(k in 1:nrow(pars_mat)){
    pdf_k[,k] = pars_mat[k,1] * dist_pdf(x, dist, pars_mat[k, -1], pmf_func)
  }

  ### summing up to get the mixture
  py <- rowSums(pdf_k, na.rm = T)

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
    warning("Some of these modes are flat.")
  }
  
  output = rep(NA_real_, length(x))
  
  if (type == "unique") {
    output[1:length(loc_modes)] = x_decrease
  }
  
  if (type == "all") {
    output[1:length(loc_modes)] = loc_modes
  }
  
  mode = list()
  mode$mode_estimates = output
  mode$dist = dist
  mode$parameters = pars
  mode$pdf_func = pmf_func
  mode$data = data
  
  class(mode) = "Mode"

  return(mode)
}