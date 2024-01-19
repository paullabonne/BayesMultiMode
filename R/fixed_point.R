#' Modal fixed-point algorithm
#' 
#' Algorithm for estimating modes in mixture of Normal distributions from Carreira-Perpinan (2000).
#' 
#' @param mixture An object of class Mixture.
#' @param tol_x Tolerance parameter for distance in-between modes; default is 1e-6; if two modes are closer than \code{tol_x}, only the first estimated mode is kept.
#' @param tol_conv Tolerance parameter for convergence of the algorithm; default is 1e-8.
#' 
#' @return Vector of estimated modes.
#' 
#' @details
#' 
#' This algorithm returns the local maxima of the mixture
#' \deqn{p(x) = \sum_{k=1}^{K}\pi_k p_k(x),}
#' where \eqn{p_k} comes from the Normal family.
#' Following Carreira-perpinan (2000), a mode \eqn{x} is found by iterating the two steps:
#' \deqn{(i) \quad p(k|x^{(n)}) = \frac{\pi_k p_k(x^{(n)})}{p(x^{(n)})},}
#' \deqn{(ii) \quad x^{(n+1)} = f(x^{(n)}),}
#' with
#' \deqn{f(x) = (\sum_k p(k|x) \sigma_k)^{-1}\sum_k p(k|x) \sigma_k \mu_k,}
#' until convergence, that is, until \eqn{abs(x^{(n+1)}-x^{(n)})< \text{tol}_\text{conv}},
#' where \eqn{\text{tol}_\text{conv}} is an argument with default value \eqn{1e-8}.
#' Following Carreira-perpinan (2000), the algorithm is started at each component location.
#' Separately, it is necessary to identify identical modes which diverge only up to
#' a small value. By default modes which are closer
#' than \eqn{sd(y)/10} are merged; this tolerance value can be controlled with the argument
#' \code{tol_x}.
#' 
#' @references
#' \insertRef{carreira-perpinan_mode-finding_2000}{BayesMultiMode}
#' 
#' @importFrom stats dnorm sd
#' @importFrom graphics abline curve
#' @importFrom assertthat assert_that
#' @importFrom stringr str_remove
#' 
#' @examples
#' mu = c(0,5)
#' sigma = c(1,2)
#' p = c(0.5,0.5)
#'
#' params = c(eta = p, mu = mu, sigma = sigma)
#' pars_names = c("eta", "mu", "sigma")
#' mix = new_Mixture(params, pars_names = pars_names)
#' modes = fixed_point(mix)
#' 
#' @export

fixed_point <- function(mixture, tol_x = 1e-6, tol_conv = 1e-8) {
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names

  ## input checks
  assert_that(length(tol_x)==1 & tol_x > 0, msg = "tol_x should be a positive scalar")
  ##
  
  modes = rep(NA_real_,length(pars)/3)
  pars = pars[!is.na(pars)]
  p = pars[grep("eta", names(pars))]
  mu = pars[grep("mu", names(pars))]
  sigma = pars[grep("sigma", names(pars))]
  
  iter = 0

  for (i in 1:length(mu)) {
    
    x = mu[i]
    delta = 1
    
    while (delta > tol_conv) {
      iter = iter + 1
      x1 = f_fp(x, p, mu, sigma)
      delta = abs(x - x1)
      x = x1
    }

    ## check that the mode is not too close to other modes
    if(any(!is.na(modes))){
      diff = abs(x-modes)
      diff = diff[!is.na(diff)]
      if (!any(diff<tol_x)) {
        modes[i] = x 
      } 
    } else {
      modes[i] = x 
    }
  }
  
  # if (show_plot) {
  #   curve(normal_mix(x, p, mu, sigma), from = min(data), to =  max(data), xlab = "", ylab = "")
  #   for (x in modes) {
  #     abline(v = x) 
  #   } 
  # }
  
  mode = list()
  mode$mode_estimates = modes
  mode$dist = mixture$dist
  mode$parameters = pars
  
  class(mode) = "Mode"
  return(mode)
}

#' @keywords internal
f_fp <- function(x, p, mu, sigma) {
  pmx = dnorm(x, mu, sigma) * p
  pmx = pmx/sum(pmx)
  
  if (any(is.na(pmx))) {
    # x yields a density of zero
    pmx = 1/length(pmx)
  }
  
  f = 1/sum(pmx/sigma^2) * sum(pmx/sigma^2*mu)
  
  return(f)
}