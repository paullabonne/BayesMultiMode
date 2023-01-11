#' Mode-finding EM algorithm (MEM)
#' 
#' From Li, Jia, Surajit Ray, and Bruce G. Lindsay.
#' "A nonparametric statistical approach to clustering via mode identification."
#' Journal of Machine Learning Research 8 (2007): 1687-1723.
#' 
#' @param mcmc a vector of estimated mixture parameters.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "student" and "skew_normal". 
#' @param tol_p Tolerance for small components. Default is 1e-3. All components with mixture weights lower than tol_p are dropped.
#' @param tol_x Tolerance for distance in-between modes. Default is sd(y)/10. If two modes are closer than tol_x, only the first estimated mode is kept.
#' @param show_plot Show the data and estimated modes.
#' 
#' @return A vector estimated modes.
#' 
#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom graphics abline curve
#' @importFrom stats optim
#' 
#' 
#' @export

MEM <- function(mcmc, dist, y, tol_p = 1e-3, tol_x, show_plot = FALSE) {
  p = mcmc[grep("theta", names(mcmc))]
  est_mode = rep(NA, length(p))
  
  mu = mcmc[grep("mu", names(mcmc))][p > tol_p]
  sigma = mcmc[grep("sigma", names(mcmc))][p > tol_p]
  
  if (dist == "student") {
    nu = mcmc[grep("nu", names(mcmc))][p > tol_p]
  } else {
    nu = NULL
  }
  
  if (dist == "skew_normal") {
    xi = mcmc[grep("xi", names(mcmc))][p > tol_p]
  } else {
    xi = NULL
  }
  
  p = p[p > tol_p]
  nK = length(p)
  post_prob = rep(NA, nK)
  
  for (j in 1:nK) {
    fnscale = -1
    
    x = mu[j]
    
    delta = 1
    
    while (delta > 1e-8) {
      # E-step
      if (dist == "skew_normal"){
        f = skew_norm_mix(x, p, mu, sigma, xi)
      }
      if (dist == "student"){
        f = student_mix(x, p, mu, sigma, nu)
      }
      
      for (k in 1:nK){
        if (dist == "skew_normal"){
          post_prob[k] = p[k] * dsn(x, xi = mu[k], omega = sigma[k], alpha = xi[k])/f 
        }
        if (dist == "student"){
          post_prob[k] = p[k] * dst(x, xi = mu[k], omega = sigma[k], nu = nu[k])/f 
        }
      }
      
      # M-step
      Min = optim(par = x, Q_func, method = "L-BFGS-B",
                  dist = dist, 
                  post_prob = post_prob,
                  mu = mu, sigma = sigma, xi = xi,
                  nu = nu,
                  control = list(fnscale = fnscale))
      
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
    
    if (x <= max(y) & x >= min(y) & not_duplicate){ #!(!(x <= mcmc["max_y"] & x >= mcmc["min_y"]) & p[j]<1e-3)
      est_mode[j] = x
    }
  }
  
  if (show_plot) {
    if (dist == "skew_normal"){
      curve(skew_norm_mix(x, p, mu, sigma, xi), from = min(y), to =  max(y))
    }
    if (dist == "student"){
      curve(student_mix(x, p, mu, sigma, nu), from = min(y), to =  max(y))
    }
    for (x in est_mode) {
      abline(v = x) 
    } 
  }
  
  return(est_mode)
}

#' @keywords internal
Q_func = function(x, dist, post_prob, mu, sigma, xi, nu, min_max = 1){
  if (dist == "student") {
    Q = sum(post_prob * log(dst_vec(x, xi = mu, omega = sigma, nu = nu)))
  }
  
  if (dist == "skew_normal") {
    Q = sum(post_prob * log(dsn(x, xi = mu, omega = sigma, alpha = xi)))
  }
  
  if(is.na(Q)|!is.finite(Q)){
    Q = -1e6
  }
  
  if(min_max==-1){
    Q = -Q
  }
  
  return(Q)
}