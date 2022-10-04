#' Density function of a mixture of shifted poisson.
#' 
#' Bayesian estimation of a mixture of shifted poisson distributions using a Sparse Finite Mixture MCMC algorithm.
#' @param x (an integer) Observation at which the density is evaluated.
#' @param p (a vector) Mixture weights.
#' @param lambda (a vector) Lambda parameter for each component.
#' @param kappa (a vector of integers) Kappa parameter for each component.
#' @returns 
#' Returns the density evaluated at the observation.
#' 
#' @examples
#' 
#' #two-component mixture
#' p = c(0.1,0.5,0.4)
#' lambda = c(1,2,3)
#' kappa= c(0,5,1)
#' dspmix(4,p,lambda,kappa)
#' dspmix(2,p,lambda,kappa)
#' dspmix(10,p,lambda,kappa)
#' 
#' @export
#'

dspmix <- function(x, p,lambda,kappa){
  mixture = 0
  for (i in 1:length(kappa)){
    mixture = mixture+p[i]*dpois(x-kappa[i], lambda[i], log = FALSE)
  }
  return(mixture)
}