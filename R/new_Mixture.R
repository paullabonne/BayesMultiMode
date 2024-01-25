#' Creating a S3 object of class \code{Mixture}
#' 
#' Function for creating an object of class \code{Mixture} which can subsequently be used as argument in [MEM()] and [fixed_point()].
#' 
#' @param pars Named vector of mixture parameters.
#' @param dist String indicating the distribution of the mixture components; default is "NA".
#' Currently supports "normal" and "skew_normal"; not needed if pdf_func is provided.
#' @param pdf_func Pdf of the mixture components; default is null.
#' @param dist_type Either "continuous" or "discrete".
#' @param data (optional) Data used for estimation; default is NULL.
#' 
#' @returns
#' A list of class \code{Mixture} containing:
#' \itemize{
#'  \item{pars}{ - Same as argument.}
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
#' dist = "skew_normal"
#' 
#' mix = new_Mixture(params, dist = dist)
#' 
#' # summary(mix)
#' # plot(mix)
#' 
#' # Example with an arbitrary distribution ===================================
#' mu = c(0,6)
#' omega = c(1,2)
#' xi = c(0,0)
#' nu = c(3,100)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = mu, sigma = omega, xi = xi, nu = nu)
##
#' 
#' pdf_func <- function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' 
#' mix = new_Mixture(params, pdf_func = pdf_func, dist_type = "continuous")
#' 
#' # summary(mix)
#' # plot(mix, from = -4, to = 4)
#' 
#' @export

new_Mixture <- function(pars,
                        dist = NA_character_,
                        pdf_func = NULL,
                        dist_type = NA_character_,
                        data = NULL) {
  ## input checks
  assert_that(is.string(dist))
  assert_that(is.string(dist_type))
  assert_that(is.vector(pars))
  assert_that(!is.null(names(pars)),
              msg = "element of pars should have names")
  assert_that(!(is.na(dist) & is.null(pdf_func)),
              msg = "one of dist or pdf_func must be specified")
  
  pars_names = unique(str_extract(names(pars), "[a-z]+"))

  assert_that("eta" %in% pars_names,
              msg = "pars should include a parameter named eta representing mixture proportions.")
  
  list_func = test_and_export(pars, pdf_func, dist, pars_names, dist_type, par_type = "pars")
  
  Mixture = list(pars = pars,
                 pars_names = pars_names,
                 dist_type = list_func$dist_type,
                 dist = dist,
                 pdf_func = list_func$pdf_func,
                 data = data,
                 nb_var = length(pars_names) - 1, #minus the shares
                 K = length(pars)/length(pars_names))
  
  class(Mixture) <- "Mixture"
  
  return(Mixture)
}