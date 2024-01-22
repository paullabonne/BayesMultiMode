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

  if (!is.na(dist)) {
    if (dist == "normal") {
      assert_that(sum(c("eta", "mu", "sigma") %in% pars_names)==3,
                  msg = "missing parameter in pars; variables should be theta, mu and sigma when using dist = normal")
    }
    
    if (dist == "skew_normal") {
      assert_that(sum(c("eta", "xi", "omega", "alpha") %in% pars_names)==4,
                  msg = "variables should be theta, xi, omega and alpha when using dist = skew_normal") 
    }
    
    if (dist == "poisson"){
      assert_that(sum(c("eta", "lambda") %in% pars_names)==2,
                  msg = "variables should be theta and lambda when using dist = poisson")
    }
    
    if (dist == "shifted_poisson"){
      assert_that(sum(c("eta", "lambda", "kappa") %in% pars_names)==3,
                  msg = "variables should be theta, lambda and kappa when using dist = shifted_poisson")
    }
    
    if (dist %in% c("normal", "skew_normal")) {
      dist_type = "continuous"
    } else if (dist %in% c("poisson", "shifted_poisson")) {
      dist_type = "discrete"
    } else {
      stop("dist must be one of the following : normal, skew_normal, poisson or shifted_poisson")
    }
  }
  
  # that pdf_func can be computed when provided
  if(!is.null(pdf_func)) {
    assert_that(!is.na(pdf_func(1, vec_to_mat(pars, pars_names)[1,-1])),
                msg = "running pdf_func with pars provided returns NA") 
    assert_that(!is.na(dist_type),
                msg = "dist_type must be provided when argument pdf_func is used") 
  }
  
  Mixture = list(pars = pars,
                 pars_names = pars_names,
                 dist_type = dist_type,
                 dist = dist,
                 pdf_func = pdf_func,
                 data = data,
                 nb_var = length(pars_names) - 1, #minus the shares
                 K = length(pars)/length(pars_names))
  
  class(Mixture) <- "Mixture"
  
  return(Mixture)
}