#' Modal Expectation-Maximization algorithm (MEM)
#' 
#' Algorithm from Li and Lindsay (2007) to find modes in mixture of continuous distributions.
#' 
#' @param mixture An object of class Mixture.
#' @param tol_x Tolerance parameter for distance in-between modes; default is 1e-6; if two modes are closer than \code{tol_x}, only the first estimated mode is kept.
#' @param tol_conv Tolerance parameter for convergence of the algorithm; default is 1e-8.
#' 
#' @return An object of class Mode.
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
#' dist = "skew_normal"
#' 
#' mix = new_Mixture(params, dist = dist)
#' modes = MEM(mix)
#' 
#' # Example with an arbitrary distribution ===================================
#' xi = c(0,6)
#' omega = c(1,2)
#' alpha = c(0,0)
#' nu = c(3,100)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = xi, sigma = omega, xi = alpha, nu = nu)
#' 
#' pdf_func <- function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' mix = new_Mixture(params, pdf_func = pdf_func, dist_type = "continuous")
#' modes = MEM(mix)
#' 
#' @export

MEM <- function(mixture, tol_x = 1e-6, tol_conv = 1e-8) {
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  dist = mixture$dist
  pdf_func = mixture$pdf_func
  
  ## input checks
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##

  pars_mat <- vec_to_mat(pars, pars_names)

  est_mode = rep(NA, nrow(pars_mat))

  nK = nrow(pars_mat)
  post_prob = rep(NA, nK)

  # remove empty components (a feature of some MCMC methods)
  pars_mat = na.omit(pars_mat)
  
  for (j in 1:nK) {
    x = pars_mat[j,2]
    
    delta = 1
    
    while (delta > 1e-8) {
      # E-step
      f_mix = dist_mixture(x, dist, pars_mat, pdf_func)
   
      for (k in 1:nK){ 
        post_prob[k] = pars_mat[k, 1] * dist_pdf(x, dist, pars_mat[k, -1], pdf_func)/f_mix 
      }
     
      # M-step
      Min = optim(par = x, Q_func, method = "L-BFGS-B",
                  dist = dist, 
                  post_prob = post_prob,
                  pars = pars_mat[, -1],
                  pdf_func = pdf_func,
                  control = list(fnscale = -1))
      
      x1 = Min$par
      
      # check convergence and increment
      delta = abs(x - x1)
      x = x1
    }

    ## check that the mode is not too close to other modes
    ## check that the mode is not too close to other modes
    if(any(!is.na(est_mode))){
      diff = abs(x-est_mode)
      diff = diff[!is.na(diff)]
      if (!any(diff<tol_x)) {
        est_mode[j] = x 
      } 
    } else {
      est_mode[j] = x 
    }
  }
  
  mode = list()
  mode$mode_estimates = est_mode[!is.na(est_mode)]
  mode$dist = dist
  mode$pars = pars
  mode$pdf_func = pdf_func
  mode$dist_type = "continuous"
  mode$algo = "Modal Expectation-Maximization (MEM)"
  mode$K = mixture$K
  mode$nb_var = mixture$nb_var
  
  class(mode) = "Mode"
  
  return(mode)
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
    stop("Q function is not finite")
  }
  
  return(Q)
}

#' @keywords internal
vec_to_mat <- function(pars, pars_names) {
  pars_mat = c()
  for (i in 1:length(pars_names)) {
    pars_mat = cbind(pars_mat, pars[grep(pars_names[i], names(pars))])
  }
  colnames(pars_mat) <- pars_names
  
  return(pars_mat)
}