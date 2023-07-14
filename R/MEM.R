#' Modal EM algorithm (MEM)
#' 
#' Algorithm from Li and Lindsay (2007) to find modes in mixture of continuous distributions.
#' 
#' @param mcmc Vector of estimated mixture parameters.
#' @param data Vector of observations used for estimating the mixture.
#' @param pars_names Names of the mixture parameters; the first element of 
#' this vector should be the name of the mixture proportions. If you have used 
#' the skew normal of Azzalini, then the second element should correspond to the location,
#' the third to the scale and the fourth to the shape.
#' @param dist String indicating the distribution of the mixture components; default is "NA".
#' Currently supports "normal" and "skew_normal"; not needed if pdf_func is provided.
#' @param pdf_func Pdf of the mixture components associated with the mcmc draws
#' (if mcmc estimation has not been carried out with BayesMultiMode); default is null.
#' @param tol_x Tolerance parameter for distance in-between modes; default is sd(data)/10; if two modes are closer than \code{tol_x}, only the first estimated mode is kept.
#' @param tol_conv Tolerance parameter for convergence of the algorithm; default is 1e-8.
#' @param show_plot If true show the data and estimated modes; default is false
#' 
#' @return Vector of estimated modes.
#' 
#' @details
#' This algorithm returns the local maxima of the mixture
#' \deqn{p(x) = \sum_{k=1}^{K}\pi_k p_k(x),}
#' where \eqn{p_k} is a density function.
#' Following Li and Lindsay (2007), a mode \eqn{x} is found by iterating the two steps:
#' \deqn{(i) \quad p(k|x^{(n)}) = \frac{\pi_k p_k(x^{(n)})}{p(x^{(n)})},}
#' \deqn{(ii) \quad x^{(n+1)} = \text{argmax}_x  \sum_k p(k|x) \text{log} p_k(x^{(n)}),}
#' until convergence, that is, until \eqn{abs(x^{(n+1)}-x^{(n)})< \text{tol}_\text{conv}},
#' where \eqn{\text{tol}_\text{conv}} is an argument with default value \eqn{1e-8}.
#' The algorithm is started at each component location.
#' Separately, it is necessary to identify identical modes which diverge only up to
#' a small value. By default modes which are closer
#' than \eqn{sd(y)/10} are merged; this tolerance value can be controlled with the argument
#' \code{tol_x}.
#' 
#' While it is also possible to use the MEM algorithm for Normal mixtures, 
#' this is not recommended because the algorithm is less efficient than the
#' fixed-point method in this particular case.
#' 
#' @references
#' \insertRef{li_nonparametric_2007}{BayesMultiMode}\cr\cr
#' 
#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom graphics abline curve
#' @importFrom stats optim na.omit var
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
#' 
#' @examples
#' 
#' # Example with the skew normal =============================================
#' xi = c(0,6)
#' omega = c(1,2)
#' alpha = c(0,0)
#' p = c(0.8,0.2)
#' params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
#' pars_names = c("eta", "xi", "omega", "alpha")
#' dist = "skew_normal"
#' 
#' data = c(sn::rsn(p[1]*100, xi[1], omega[1], alpha[1]),
#'          sn::rsn(p[2]*100, xi[2], omega[2], alpha[2]))
#' 
#' modes = MEM(params, data = data, pars_names = pars_names, dist = dist)
#' 
#' # Example with an arbitrary distribution ===================================
#' xi = c(0,6)
#' omega = c(1,2)
#' alpha = c(0,0)
#' nu = c(3,100)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = xi, sigma = omega, xi = alpha, nu = nu)
#' pars_names = c("eta", "mu", "sigma", "xi", "nu")
#' 
#' pdf_func <- function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' data = c(sn::rst(p[1]*100, xi[1], omega[1], alpha = alpha[1], nu = nu[1]),
#'          sn::rst(p[2]*100, xi[2], omega[2], alpha = alpha[2], nu = nu[2]))
#' 
#' modes = MEM(params, pars_names = pars_names, data = data, pdf_func = pdf_func)
#' 
#' @export

MEM <- function(mcmc, data, pars_names, dist = "NA", pdf_func = NULL, tol_x = sd(data)/10, tol_conv = 1e-8, show_plot = FALSE) {
  ## input checks
  fail = "inputs to the Mode-finding EM algorithm are corrupted"
  assert_that(is.vector(mcmc) & length(mcmc) >= 3,
              msg = "mcmc should be a vector of length >= 3")
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.logical(show_plot), msg = "show_plot should be TRUE or FALSE")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##
  
  ##
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
              msg = "the name of the parameters provided by pars_names and those of the mcmc vector do not match")
  
  if (dist %in% c("skew_normal")) {
    assert_that(length(pars_names) == 4,
                msg = "the number of elements in pars_names does not match with dist") 
  }
  ##
  
  pars = c()
  for (i in 1:length(pars_names)) {
    pars = cbind(pars, mcmc[grep(pars_names[i], names(mcmc))])
  }
  
  colnames(pars) <- pars_names

  est_mode = rep(NA, nrow(pars))

  nK = nrow(pars)
  post_prob = rep(NA, nK)

  # remove empty components (a feature of some MCMC methods)
  pars = na.omit(pars)
  
  for (j in 1:nK) {
    x = pars[j,2]
    
    delta = 1
    
    while (delta > 1e-8) {
      # E-step
      f_mix = dist_mixture(x, dist, pars, pdf_func)
      
      for (k in 1:nK){ 
        post_prob[k] = pars[k, 1] * dist_pdf(x, dist, pars[k, -1], pdf_func)/f_mix 
      }
     
      # M-step
      Min = optim(par = x, Q_func, method = "L-BFGS-B",
                  dist = dist, 
                  post_prob = post_prob,
                  pars = pars[, -1],
                  pdf_func = pdf_func,
                  control = list(fnscale = -1))
      
      x1 = Min$par
      
      # check convergence and increment
      delta = abs(x - x1)
      x = x1
    }
 
    ## check that the mode is not too close to other modes
    not_duplicate = TRUE
    
    if (any(is.finite(est_mode))) {
      diff = abs(x-est_mode)
      diff = diff[!is.na(diff)]
      if (any(diff<tol_x)) {
        not_duplicate = FALSE
      }
    }
    
    if (x <= max(data) & x >= min(data) & not_duplicate){
      est_mode[j] = x
    }
  }
  
  if (show_plot) {
    curve(dist_mixture(x, dist, pars, pdf_func), from = min(data), to =  max(data))
    for (x in est_mode) {
      abline(v = x) 
    } 
  }
  
  return(est_mode)
}

#' @keywords internal
Q_func = function(x, dist, post_prob, pars, pdf_func){
  
  pdf = rep(NA, nrow(pars))
  
  for (i in 1:nrow(pars)) {
    pdf[i] = dist_pdf(x, dist, pars[i,], pdf_func = pdf_func)
  } 
  
  pdf[pdf==0] = 1e-10 #otherwise the log operation below can return infs
  
  Q = sum(post_prob * log(pdf))
  
  if(is.na(Q)|!is.finite(Q)){
    Q = -1e6
  }
  
  return(Q)
}