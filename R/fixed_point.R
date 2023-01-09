#' Fixed-point algorithm for finding the modes of a gaussian mixture
#' 
#' Fixed-point algorithm to find the modes of a gaussian
#' mixture from Carreira-Perpinan (2000), section 4 equation (10)
#' https://doi.org/10.1109/34.888716
#' 
#' @importFrom stats dnorm
#' @importFrom graphics abline curve
#' @export

fixed_point <- function(params, tol_p = 1e-4, tol_x, show_plot = FALSE) {
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
    
    if (x <= params["max_y"] & x >= params["min_y"] & not_duplicate){
      modes[i] = x 
    }
  }
  
  if (show_plot) {
    curve(Gaussian_mixture(x, p, mu, sigma), from = params["min_y"], to =  params["max_y"])
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