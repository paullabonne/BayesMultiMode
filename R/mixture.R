#' Creating a S3 object of class `mixture`
#' 
#' Creates an object of class `mixture` which can subsequently be used as argument in [mix_mode()] for mode estimation.
#' 
#' @param pars Named vector of mixture parameters.
#' @param dist Distribution family of the mixture components supported by
#' the package (i.e. `"normal"`, `"student"`, `"skew_normal"` or `"shifted_poisson"`).
#' If left unspecified, `pdf_func` is required.
#' @param pdf_func (function) Pdf or pmf of the mixture components;
#' this input is used only if `dist` is left unspecified.
#' pdf_func should have two arguments : (i) the observation where the pdf is evaluated;
#' (ii) a named vector representing the function parameters. For instance a normal pdf would take the form:
#' `pdf_func <- function(x, par) dnorm(x, par['mu'], par['sigma'])`.
#' The names of `par` should correspond to variables in `pars`, e.g. `"mu1"`, `"mu2"` etc... 
#' @param dist_type Type of the distribution, either `"continuous"` or `"discrete"`.
#' @param range (optional for continuous mixtures) upper and lower limit of the range where the mixture should be evaluated.
#' @param loc (for continuous mixtures other than Normal mixtures) String indicating the location parameter
#' of the distribution; the latter is used to initialise the MEM algorithm.
#' 
#' @returns
#' A list of class `mixture` containing:
#'  \item{pars}{Same as argument.}
#'  \item{pars_names}{Names of the parameters of the components' distribution.}
#'  \item{dist}{Same as argument.}
#'  \item{pdf_func}{Pdf (or pmf) of the mixture components.}
#'  \item{dist_type}{Same as argument.}
#'  \item{loc}{Type of the distribution, either `"continuous"` or `"discrete"`.}
#'  \item{nb_var}{Number of parameters in the mixture distribution.}
#'  \item{K}{Number of mixture components.}
#'  \item{range}{Same as argument.}
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
#' mix = mixture(params, dist = dist)
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
#' mix = mixture(params, pdf_func = pdf_func,
#' dist_type = "continuous", loc = "mu")
#' 
#' # summary(mix)
#' # plot(mix, from = -4, to = 4)
#' 
#' @export

mixture <- function(pars,
                        dist = NA_character_,
                        pdf_func = NULL,
                        dist_type = NA_character_,
                        range = NULL,
                        loc = NA_character_) {
  ## input checks
  assert_that(is.string(dist))
  assert_that(is.string(dist_type))
  assert_that(is.vector(pars))
  assert_that(!is.null(names(pars)),
              msg = "element of pars should have names")
  assert_that(!(is.na(dist) & is.null(pdf_func)),
              msg = "one of dist or pdf_func must be specified")
  
  pars_names = unique(str_extract(names(pars), "[a-z]+"))
  
  list_func = test_and_export(pars, pdf_func, dist, pars_names, dist_type, loc)
  
  if (list_func$dist_type == "discrete") {
    assert_that(!is.null(range),
                msg = "range argument must be filled when using a discrete distribution")
    assert_that(is.vector(range) & length(range) == 2,
                msg = "range should be a vector of length 2")
    assert_that(all(is.finite(range)),
                msg = "lower and upper limits of range should be finite")
    assert_that(range[2] > range[1],
                msg = "upper limit of range not greater than lower limit")
    
    if (dist %in% c("poisson", "shifted_poisson")) {
      assert_that(all(range>=0),
                  msg = "lower limit should be greater or equal than zero when using the Poisson or shifted Poisson.")
    }
  }
  
  mixture = list(pars = pars,
                 pars_names = pars_names,
                 dist_type = list_func$dist_type,
                 dist = dist,
                 pdf_func = list_func$pdf_func,
                 range = range,
                 loc = list_func$loc,
                 nb_var = length(pars_names) - 1, #minus the shares
                 K = list_func$K)
  
  class(mixture) <- "mixture"
  
  return(mixture)
}