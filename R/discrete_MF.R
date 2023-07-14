#' Mode-finding algorithm for mixture of discrete distributions
#' 
#' Function to estimate modes in mixtures of discrete distributions following; see Cross et al. (2023).
#' 
#' @param mcmc Vector of estimated mixture parameters.
#' @param data Vector of observations used for estimating the mixture.
#' @param pars_names Names of the mixture parameters; first element should 
#' correspond to the mixture proportions. 
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "poisson" and "shifted_poisson"; default is "NA"; only
#' use this argument if you have used Poisson and shifted Poisson distributions
#' identical to the one used in the package.
#' @param pmf_func Pmf of the mixture components associated with the mcmc draws.
#' (if mcmc estimation has not been carried out with \pkg{BayesMultiMode}); default is null.
#' @param type Type of modes, either unique or all (the latter includes flat modes); default is "all".
#' @param show_plot If true show the data and estimated modes; default is false.
#' 
#' @return Vector of estimated modes.
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
#' pars_names = c("eta", "lambda")
#' dist = "poisson"
#' 
#' data = c(rpois(p[1]*1e3, lambda[1]),
#'          rpois(p[2]*1e3, lambda[2]))
#' 
#' modes = discrete_MF(params, data = data, pars_names = pars_names, dist = dist)
#' 
#' # Example with an arbitrary distribution ===================================
#' mu = c(20,5)
#' size = c(20,0.5)
#' p = c(0.5,0.5)
#' params = c(eta = p, mu = mu, size = size)
#' pars_names = c("eta", "mu", "size")
#' 
#' data = c(rnbinom(p[1]*1e3, mu = mu[1], size = size[1]),
#'          rnbinom(p[2]*1e3, mu = mu[2], size = size[2]))
#' 
#' pmf_func <- function(x, pars) {
#'   dnbinom(x, mu = pars["mu"], size = pars["size"])
#' }
#' 
#' modes = discrete_MF(params, data = data, pars_names = pars_names, pmf_func = pmf_func)
#' 
#' @export

discrete_MF <- function(mcmc, data, pars_names, dist = "NA",
                        pmf_func = NULL, type = "all", show_plot = FALSE){
  
  ## input checks
  assert_that(is.vector(mcmc),
              msg = "mcmc should be a vector")
  assert_that(is.string(dist),
              msg = "dist should be a string")
  
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(!any(is.na(data)) & !any(is.infinite(data)),
              msg = "y should not include missing or infinite values")
  assert_that(type %in% c("unique", "all"),
              msg = "type should be either 'unique' or 'all' ")
  assert_that(is.logical(show_plot), msg = "show_plot should be either TRUE or FALSE")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##
  
  ##
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  x = min(data):max(data)
  
  assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
              msg = "the name of the parameters provided by pars_names and those of the mcmc vector do not match")
  
  if (dist %in% c("shifted_poisson")) {
    assert_that(length(pars_names) == 3,
                msg = "the number of elements in pars_names does not match with dist")
  }
  if (dist %in% c("poisson")) {
    assert_that(length(pars_names) == 2,
                msg = "the number of elements in pars_names does not match with dist")
  }
  ##
  
  pars = c()
  for (i in 1:length(pars_names)) {
    pars = cbind(pars, mcmc[grep(pars_names[i], names(mcmc))])
  }
  
  colnames(pars) <- pars_names
  ##
  
  Khat = nrow(pars)
  
  ### Getting individual component densities
  
  pdf_k = matrix(0, nrow=length(x), ncol=Khat) 
  for(k in 1:nrow(pars)){
    pdf_k[,k] = pars[k,1] * dist_pdf(x, dist, pars[k, -1], pmf_func)
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
    cat("Warning : Some of these modes are flat.")
  }
  
  if (show_plot) {
    plot(x, py, type = "l", xlab = "", ylab = "")
    for (mode_i in loc_modes) {
      abline(v = mode_i)
    }
  }
  
  output = rep(NA, length(x))
  
  if (type == "unique") {
    output[1:length(loc_modes)] = x_decrease
  }
  
  if (type == "all") {
    output[1:length(loc_modes)] = loc_modes
  }

  return(output)
}