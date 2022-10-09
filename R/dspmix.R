#' Density of a mixture of shifted Poisson distributions
#' 
#' Density function of a mixture of shifted Poisson distributions.
#' @param x (an integer) Observation at which the density is evaluated.
#' @param p (a vector) Mixture weights.
#' @param lambda (a vector) Lambda parameter for each component.
#' @param kappa (a vector of integers) Kappa parameter for each component.
#' @returns 
#' Returns the density evaluated at the observation.
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{basturk_bayes_2021}{BayesMultiMode}
#' @examples
#' 
#' # a three-component mixture
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