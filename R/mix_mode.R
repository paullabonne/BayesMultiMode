#' Find modes of mixtures
#' 
#' This function estimate modes in univariate mixture distributions.
#' The fixed-point algorithm of Carreira-Perpinan (2000) is used for Gaussian mixtures.
#' The Modal EM algorithm of Li et al. (2007) is used for other continuous mixtures.
#' A basic algorithm is used for discrete mixtures (see Cross et al. 2023).
#' 
#' @param mixture An object of class \code{Mixture} generated with \code{new_Mixture()}.
#' @param tol_x (for continuous mixtures) Tolerance parameter for distance in-between modes; default is 1e-6; if two modes are closer than \code{tol_x}, only the first estimated mode is kept.
#' @param tol_conv (for continuous mixtures) Tolerance parameter for convergence of the algorithm; default is 1e-8.
#' @param type (for discrete mixtures) Type of modes, either unique or all (the latter includes flat modes); default is "all".
#' 
#' @return An object of class Mode.
#' 
#' @references
#' \insertRef{cross_2023}{BayesMultiMode}
#' 
#' @details
#' 
#' This function finds modes in a univariate mixture defined as:
#' \deqn{p(.) = \sum_{k=1}^{K}\pi_k p_k(.),}
#' where \eqn{p_k} is a density or probability mass function.
#' 
#' \cr\cr
#' Fixed-point algorithm\cr
#' Following Carreira-perpinan (2000), a mode \eqn{x} is found by iterating the two steps:
#' \deqn{(i) \quad p(k|x^{(n)}) = \frac{\pi_k p_k(x^{(n)})}{p(x^{(n)})},}
#' \deqn{(ii) \quad x^{(n+1)} = f(x^{(n)}),}
#' with
#' \deqn{f(x) = (\sum_k p(k|x) \sigma_k)^{-1}\sum_k p(k|x) \sigma_k \mu_k,}
#' until convergence, that is, until \eqn{abs(x^{(n+1)}-x^{(n)})< \text{tol}_\text{conv}},
#' where \eqn{\text{tol}_\text{conv}} is an argument with default value \eqn{1e-8}.
#' Following Carreira-perpinan (2000), the algorithm is started at each component location.
#' Separately, it is necessary to identify identical modes which diverge only up to
#' a small value. By default modes which are closer
#' than \eqn{sd(y)/10} are merged; this tolerance value can be controlled with the argument
#' \code{tol_x}.
#' 
#' \cr\cr
#' MEM algorithm\cr
#' Following Li and Lindsay (2007), a mode \eqn{x} is found by iterating the two steps:
#' \deqn{(i) \quad p(k|x^{(n)}) = \frac{\pi_k p_k(x^{(n)})}{p(x^{(n)})},}
#' \deqn{(ii) \quad x^{(n+1)} = \text{argmax}_x  \sum_k p(k|x) \text{log} p_k(x^{(n)}),}
#' until convergence, that is, until \eqn{abs(x^{(n+1)}-x^{(n)})< \text{tol}_\text{conv}},
#' where \eqn{\text{tol}_\text{conv}} is an argument with default value \eqn{1e-8}.
#' The algorithm is started at each component location.
#' Separately, it is necessary to identify identical modes which diverge only up to
#' a small value. By default modes which are closer
#' than \eqn{sd(y)/10} are merged; this tolerance value can be controlled with the argument
#' \code{tol_x}.
#' 
#' \cr\cr
#' Discrete method\cr
#' By definition, modes must satisfy either: 
#'  \deqn{p_k(y_{m}-1) < p_k(y_{m}) > p_k(y_{m}+1)};
#'  \deqn{p_k(y_{m}-1) < p_k(y_{m}) = p_k(y_{m}+1) = \ldots = p_k(y_{m}+l-1) > p_k(y_{m}+l).}
#'  
#'  The algorithm evaluate each location point with these two conditions.
#' 
#' @references
#' \insertRef{li_nonparametric_2007}{BayesMultiMode}\cr\cr
#' \insertRef{carreira-perpinan_mode-finding_2000}{BayesMultiMode}
#' \insertRef{schaap_genome-wide_2013}{BayesMultiMode}
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
#' @importFrom stringr str_remove
#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom graphics abline curve
#' @importFrom stats dnorm sd optim na.omit var
#' 
#' @examples
#' 
#' # Example with a normal distribution ====================================
#' mu = c(0,5)
#' sigma = c(1,2)
#' p = c(0.5,0.5)
#'
#' params = c(eta = p, mu = mu, sigma = sigma)
#' mix = new_Mixture(params, dist = "normal")
#' modes = mix_mode(mix)
#' 
#' # summary(modes)
#' # plot(modes)
#' 
#' # Example with a skew normal =============================================
#' xi = c(0,6)
#' omega = c(1,2)
#' alpha = c(0,0)
#' p = c(0.8,0.2)
#' params = c(eta = p, xi = xi, omega = omega, alpha = alpha)
#' dist = "skew_normal"
#' 
#' mix = new_Mixture(params, dist = dist)
#' modes = mix_mode(mix)
#' # summary(modes)
#' # plot(modes)
#' 
#' # Example with an arbitrary continuous distribution ======================
#' xi = c(0,6)
#' omega = c(1,2)
#' alpha = c(0,0)
#' nu = c(3,100)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = xi, sigma = omega, xi = alpha, nu = nu)
#' 
#' pdf_func <- function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' mix = new_Mixture(params, pdf_func = pdf_func, dist_type = "continuous")
#' modes = mix_mode(mix)
#' 
#' # summary(modes)
#' # plot(modes, from = -4, to = 4)
#' 
#' # Example with a poisson distribution ====================================
#' lambda = c(0.1,10)
#' p = c(0.5,0.5)
#' params = c(eta = p, lambda = lambda)
#' dist = "poisson"
#' 
#' data = c(rpois(p[1]*1e3, lambda[1]),
#'          rpois(p[2]*1e3, lambda[2]))
#' 
#' mix = new_Mixture(params, data = data, dist = dist)
#' 
#' modes = mix_mode(mix)
#'
#' # summary(modes)
#' # plot(modes, from = 0, to = 20)
#' 
#' # Example with an arbitrary discrete distribution =======================
#' mu = c(20,5)
#' size = c(20,0.5)
#' p = c(0.5,0.5)
#' params = c(eta = p, mu = mu, size = size)
#' 
#' data = c(rnbinom(p[1]*1e3, mu = mu[1], size = size[1]),
#'          rnbinom(p[2]*1e3, mu = mu[2], size = size[2]))
#' 
#' pmf_func <- function(x, pars) {
#'   dnbinom(x, mu = pars["mu"], size = pars["size"])
#' }
#' 
#' mix = new_Mixture(params, data = data,
#' pdf_func = pmf_func, dist_type = "discrete")
#' modes = mix_mode(mix)
#' 
#' # summary(modes)
#' # plot(modes, from = 0, to = 50)
#' 
#' @export

mix_mode <- function(mixture, tol_x = 1e-6, tol_conv = 1e-8, type = "all") {
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  assert_that(all(c("pars", "pars_names", "dist_type",
                "dist", "pdf_func", "data", "nb_var", "K") %in% names(mixture)),
              msg = "mixture object is missing arguments.") 
  assert_that(length(tol_x)==1 & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(length(tol_conv)==1 & tol_conv > 0, msg = "tol_conv should be a positive scalar")
  
  pars = mixture$pars
  pars_names = mixture$pars_names
  dist = mixture$dist
  dist_type = mixture$dist_type
  ## input checks
  assert_that(length(tol_x)==1 & tol_x > 0, msg = "tol_x should be a positive scalar")
  ##
  
  if (dist_type == "continuous") {
    if (!is.na(dist) && dist == "normal") {
      modes = fixed_point(mixture, tol_x, tol_conv)
    } else {
      modes = MEM(mixture, tol_x, tol_conv)
    }
  }
  
  if (dist_type == "discrete") {
      modes = discrete_MF(mixture)
  }
  
  return(modes)
}

### internal functions
#' @keywords internal
fixed_point <- function(mixture, tol_x = 1e-6, tol_conv = 1e-8) {
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  
  ## input checks
  assert_that(mixture$dist == "normal", msg = "fixed_point only works for normal mixtures (dist = normal).")
  ##
  
  modes = rep(NA_real_,length(pars)/3)
  pars = pars[!is.na(pars)]
  p = pars[grep("eta", names(pars))]
  mu = pars[grep("mu", names(pars))]
  sigma = pars[grep("sigma", names(pars))]
  
  iter = 0
  
  for (i in 1:length(mu)) {
    
    x = mu[i]
    delta = 1
    
    while (delta > tol_conv) {
      iter = iter + 1
      x1 = f_fp(x, p, mu, sigma)
      delta = abs(x - x1)
      x = x1
    }
    
    ## check that the mode is not too close to other modes
    if(any(!is.na(modes))){
      diff = abs(x-modes)
      diff = diff[!is.na(diff)]
      if (!any(diff<tol_x)) {
        modes[i] = x 
      } 
    } else {
      modes[i] = x 
    }
  }
  
  mode = list()
  mode$mode_estimates = modes[!is.na(modes)]
  mode$dist = mixture$dist
  mode$pars = pars
  mode$dist_type = "continuous"
  mode$algo = "fixed-point"
  mode$K = mixture$K
  mode$nb_var = mixture$nb_var
  
  class(mode) = "Mode"
  return(mode)
}

#' @keywords internal
f_fp <- function(x, p, mu, sigma) {
  pmx = dnorm(x, mu, sigma) * p
  pmx = pmx/sum(pmx)
  
  if (any(is.na(pmx))) {
    # x yields a density of zero
    pmx = 1/length(pmx)
  }
  
  f = 1/sum(pmx/sigma^2) * sum(pmx/sigma^2*mu)
  
  return(f)
}

#' @keywords internal
MEM <- function(mixture, tol_x = 1e-6, tol_conv = 1e-8) {
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  pdf_func = mixture$pdf_func
  
  ## input checks
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.vector(pars_names) & is.character(pars_names),
              msg = "pars_names should be a character vector")
  ##
  
  pars_mat <- vec_to_mat(pars, pars_names)
  
  est_mode = rep(NA, nrow(pars_mat))
  
  nK = nrow(pars_mat)
  post_prob = rep(NA, nK)
  
  # remove empty components (a feature of some MCMC methods)
  pars_mat = na.omit(pars_mat)
  
  for (j in 1:nK) {
    x = pars_mat[j,2]
    
    delta = 1
    
    while (delta > 1e-8) {
      # E-step
      f_mix = pdf_func_mix(x, pars_mat, pdf_func)
      
      for (k in 1:nK){ 
        post_prob[k] = pars_mat[k, 1] * pdf_func(x, pars_mat[k, -1])/f_mix 
      }
      
      # M-step
      Min = optim(par = x, Q_func, method = "L-BFGS-B",
                  post_prob = post_prob,
                  pars = pars_mat[, -1],
                  pdf_func = pdf_func,
                  control = list(fnscale = -1))
      
      x1 = Min$par
      
      # check convergence and increment
      delta = abs(x - x1)
      x = x1
    }
    
    ## check that the mode is not too close to other modes
    ## check that the mode is not too close to other modes
    if(any(!is.na(est_mode))){
      diff = abs(x-est_mode)
      diff = diff[!is.na(diff)]
      if (!any(diff<tol_x)) {
        est_mode[j] = x 
      } 
    } else {
      est_mode[j] = x 
    }
  }
  
  mode = list()
  mode$mode_estimates = est_mode[!is.na(est_mode)]
  mode$dist = mixture$dist
  mode$pars = pars
  mode$pdf_func = pdf_func
  mode$dist_type = "continuous"
  mode$algo = "Modal Expectation-Maximization (MEM)"
  mode$K = mixture$K
  mode$nb_var = mixture$nb_var
  
  class(mode) = "Mode"
  
  return(mode)
}

#' @keywords internal
Q_func = function(x, post_prob, pars, pdf_func){
  
  pdf = rep(NA, nrow(pars))
  
  for (i in 1:nrow(pars)) {
    pdf[i] = pdf_func(x, pars[i,])
  } 
  
  pdf[pdf==0] = 1e-10 #otherwise the log operation below can return infs
  
  Q = sum(post_prob * log(pdf))
  
  if(is.na(Q)|!is.finite(Q)){
    stop("Q function is not finite")
  }
  
  return(Q)
}

#' @keywords internal
discrete_MF <- function(mixture, type = "all"){
  assert_that(inherits(mixture, "Mixture"), msg = "mixture should be an object of class Mixture")
  pars = mixture$pars
  pars_names = mixture$pars_names
  dist = mixture$dist
  pdf_func = mixture$pdf_func
  data = mixture$data
  
  ## input checks
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(!any(is.na(data)) & !any(is.infinite(data)),
              msg = "data should not include missing or infinite values")
  assert_that(type %in% c("unique", "all"),
              msg = "type must be either 'unique' or 'all' ")
  ##
  
  ##
  x = min(data):max(data)
  ##
  
  pars_mat <- vec_to_mat(pars, pars_names)
  ##
  
  Khat = nrow(pars_mat)
  
  ### Getting denisty
  py = pdf_func_mix(x, pars_mat, pdf_func)
  
  # change in the pdf
  d_py = diff(py)
  
  # where does the pdf decrease ?
  x_decrease = x[d_py<0]
  
  # Only keep the points where the pdf starts to decrease; these are modes
  d2_py = c(0, x_decrease[-1] - x_decrease[-length(x_decrease)])
  x_decrease = x_decrease[which(d2_py!=1)]
  
  # get pdf at these modes 
  pdf_modes = py[x %in% x_decrease]
  
  # get all the points at these peaks (there might be flat modes) 
  loc_modes = x[which(py %in% pdf_modes)]
  
  if (length(loc_modes) != length(x_decrease)) {
    warning("Some modes are flat.")
  }
  
  output = rep(NA_real_, length(x))
  
  if (type == "unique") {
    output[1:length(loc_modes)] = x_decrease
  }
  
  if (type == "all") {
    output[1:length(loc_modes)] = loc_modes
  }
  
  mode = list()
  mode$mode_estimates = output[!is.na(output)]
  mode$dist = dist
  mode$pars = pars
  mode$pdf_func = pdf_func
  mode$data = data
  mode$dist_type = "discrete"
  mode$py = py
  mode$algo = "discrete"
  mode$K = mixture$K
  mode$nb_var = mixture$nb_var
  
  class(mode) = "Mode"
  
  return(mode)
}