#' Fixed-point algorithm for finding the modes of a gaussian mixture.
#' See Carreira-Perpinan (2000), section 4 equation (10) https://doi.org/10.1109/34.888716.
#' 
#' @param params a vector of estimated mixture parameters.
#' @param y the vector of data used for estimating the mixtures.
#' @param tol_p Tolerance for small components. Default is 1e-3. All components with mixture weights lower than tol_p are dropped.
#' @param tol_x Tolerance for distance in-between modes. Default is sd(y)/10. If two modes are closer than tol_x, only the first estimated mode is kept.
#' @param show_plot Show the data and estimated modes.
#' 
#' @return A vector of estimated modes. 
#' 
#' @importFrom stats dnorm sd 
#' @importFrom graphics abline curve
#' @export

fixed_point <- function(params, y, tol_p = 1e-3, tol_x = sd(y)/10, show_plot = FALSE) {
  p = params[grep("theta", names(params))]
  mu = params[grep("mu", names(params))][p > tol_p]
  sigma = params[grep("sigma", names(params))][p > tol_p]
  
  modes = rep(NA,length(p))
  p = p[p > tol_p]
  
  iter = 0
  
  for (i in 1:length(mu)) {
    x = mu[i]
    delta = 1
    
    while (delta > 1e-20) {
      iter = iter + 1
      x1 = f_fp(x, p, mu, sigma)
      x = x1
      delta = abs(x - x1)
    }
    
    ## check that the mode is not too close to other modes
    not_duplicate = TRUE
    
    if (any(is.finite(modes))) {
      diff = abs(x-modes)
      diff = diff[!is.na(diff)]
      if (any(diff<tol_x)) {
        not_duplicate = FALSE
      }
    }
    
    if (x <= max(y) & x >= min(y) & not_duplicate){
      modes[i] = x 
    }
  }
  
  if (show_plot) {
    curve(normal_mix(x, p, mu, sigma), from = min(y), to =  max(y))
    for (x in modes) {
      abline(v = x) 
    } 
  }
  
  return(modes)
}

#' @keywords internal
f_fp <- function(x, p, mu, sigma) {
  pmx = dnorm(x, mu, sigma) * p
  pmx = pmx / sum(pmx)
  
  f = 1 / sum(pmx / sigma ^ 2) * sum(pmx / sigma ^ 2 * mu)
  
  return(f)
}