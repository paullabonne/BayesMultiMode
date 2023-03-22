#' Fixed-point algorithm for finding the modes of a gaussian mixture.
#' 
#' See Carreira-Perpinan (2000), section 4.
#' 
#' @param mcmc Vector of estimated mixture parameters.
#' @param data Vector of observations used for estimating the mixture.
#' @param tol_x Tolerance parameter for distance in-between modes. Default is sd(data)/10. If two modes are closer than tol_x, only the first estimated mode is kept.
#' @param show_plot If true show the data and estimated modes; default is false.
#' 
#' @return Vector of estimated modes. 
#' 
#' @references
#' \insertRef{carreira-perpinan_mode-finding_2000}{BayesMultiMode}
#' 
#' @importFrom stats dnorm sd
#' @importFrom graphics abline curve
#' @importFrom assertthat assert_that
#' @importFrom stringr str_remove
#' @export

fixed_point <- function(mcmc, data, tol_x = sd(data)/10, show_plot = F) {
  
  ## input checks
  fail = "inputs to the fixed point algorithm are corrupted"
  assert_that(is.vector(mcmc) & length(mcmc) >= 3,
              msg = paste0("mcmc should be a vector of length >= 3; ", fail))
  assert_that(is.vector(data) & length(data) > 0,
              msg = paste0("data should be a vector of length > 0; ", fail))
  assert_that(is.vector(tol_x) & tol_x > 0, msg = paste0("tol_x should be a positive scalar; ", fail))
  assert_that(is.logical(show_plot), msg = paste0("show_plot should be TRUE or FALSE; ", fail))
  
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(c("theta", "mu", "sigma") %in% names_mcmc)==3,
              msg = paste0("missing parameter in mcmc; ", fail))
  ##
  
  modes = rep(NA,length(mcmc)/3)
  mcmc = mcmc[!is.na(mcmc)]
  p = mcmc[grep("theta", names(mcmc))]
  mu = mcmc[grep("mu", names(mcmc))]
  sigma = mcmc[grep("sigma", names(mcmc))]

  assert_that(length(p) == length(mu) & length(sigma) == length(mu),
              msg = paste0("p, mu and sigma should have the same lengths", fail))
  
  iter = 0
  
  
  
  for (i in 1:length(mu)) {
    
    x = mu[i]
    delta = 1
    
    while (delta > 1e-8) {
      iter = iter + 1
      x1 = f_fp(x, p, mu, sigma)
      delta = abs(x - x1)
      x = x1
      
      if (is.na(delta)){browser()}
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
    
    if (x <= max(data) & x >= min(data) & not_duplicate){
      modes[i] = x 
    }
  }
  
  if (show_plot) {
    curve(normal_mix(x, p, mu, sigma), from = min(data), to =  max(data), xlab = "", ylab = "")
    for (x in modes) {
      abline(v = x) 
    } 
  }
  
  return(modes)
}

#' @keywords internal
f_fp <- function(x, p, mu, sigma) {
  pmx = dnorm(x, mu, sigma) * p
  pmx = pmx/sum(pmx)

  f = 1/sum(pmx/sigma^2) * sum(pmx/sigma^2*mu)
  
  if (is.na(f)){browser()}
  
  return(f)
}