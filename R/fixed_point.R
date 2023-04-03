#' Modal fixed-point algorithm
#' 
#' Algorithm for estimating modes in mixture of Normal distributions.
#' 
#' @param mcmc Vector of estimated mixture parameters
#' @param data Vector of observations used for estimating the mixture
#' @param pars_names Names of the mixture parameters; first element should 
#' correspond to the mixture proportions; second to the mean; third to the 
#' standard deviation.
#' @param tol_x Tolerance parameter for distance in-between modes; default is sd(data)/10; if two modes are closer than tol_x, only the first estimated mode is kept.
#' @param show_plot If true show the data and estimated modes; default is false
#' 
#' @return Vector of estimated modes 
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
#' data = c(rnorm(p[1]*100, mu[1], sigma[1]), rnorm(p[2]*100, mu[2], sigma[2]))
#' params = c(eta = p, mu = mu, sigma = sigma)
#' pars_names = c("eta", "mu", "sigma")
#' modes = fixed_point(params, data, pars_names)
#' 
#' @export

fixed_point <- function(mcmc, data, pars_names, tol_x = sd(data)/10, show_plot = F) {
  
  ## input checks
  assert_that(is.vector(mcmc) & length(mcmc) >= 3,
              msg = "mcmc should be a vector of length >= 3; ")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0; ")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.logical(show_plot), msg = "show_plot should be TRUE or FALSE")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(c("eta", "mu", "sigma") %in% names_mcmc)==3,
              msg = "missing parameter in mcmc; variables should be theta, mu and sigma")
  ##
  
  modes = rep(NA,length(mcmc)/3)
  mcmc = mcmc[!is.na(mcmc)]
  p = mcmc[grep(pars_names[1], names(mcmc))]
  mu = mcmc[grep(pars_names[2], names(mcmc))]
  sigma = mcmc[grep(pars_names[3], names(mcmc))]

  assert_that(length(p) == length(mu) & length(sigma) == length(mu),
              msg = "p, mu and sigma should have the same lengths")
  
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

  return(f)
}