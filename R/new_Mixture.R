#' Creating a S3 object of class \code{Mixture}
#' 
#' Function for creating an object of class \code{Mixture} which can subsequently be used as argument in [MEM()] and [fixed_point()].
#' 
#' @param pars Vector of mixture parameters.
#' @param pars_names Names of the mixture parameters; the first element of 
#' this vector should be the name of the mixture proportions. If you have used 
#' the skew normal of Azzalini, then the second element should correspond to the location,
#' the third to the scale and the fourth to the shape.
#' @param dist String indicating the distribution of the mixture components; default is "NA".
#' Currently supports "normal" and "skew_normal"; not needed if pdf_func is provided.
#' @param pdf_func Pdf of the mixture components; default is null.
#' @param dist_type Either "continuous" or "discrete".
#' 
#' @returns
#' A list of class \code{Mixture} containing:
#' \itemize{
#'  \item{pars}{ - Same as argument.}
#'  \item{pars_names}{ - Same as argument.}
#'  \item{dist}{ - Same as argument.}
#'  \item{pdf_func}{ - Same as argument.}
#'  \item{dist_type}{ - Same as argument.}
#' }
#' 
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
#' 
#' mix = new_Mixture(params, pars_names = pars_names, dist = dist)
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
#' 
#' mix = new_Mixture(params, pars_names = pars_names, pdf_func = pdf_func)
#' 
#' @export

new_Mixture <- function(pars,
                        pars_names,
                        dist_type = "NA",
                        dist = "NA",
                        pdf_func = NULL) {
  ## input checks
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.string(dist_type),
              msg = "dist_type should be a string")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##
  
  Mixture = list(pars = pars,
                 pars_names = pars_names,
                 dist_type = dist_type,
                 dist = dist,
                 pdf_func = pdf_func)
  
  class(Mixture) <- "Mixture"
  
  return(Mixture)
}