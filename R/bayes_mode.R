#' Bayesian mode inference
#' 
#' Bayesian inference on the modes in a univariate mixture estimated with MCMC methods, see \insertCite{Cross2024;textual}{BayesMultiMode}.
#' Provides posterior probabilities of the number of modes and their locations.
#' Under the hood it calls the function [mix_mode()] to find the modes in each MCMC draw.
#' 
#' @param BayesMix An object of class `bayes_mixture` generated with either [bayes_fit()] or [bayes_mixture()].
#' @param rd (for continuous mixtures) Integer indicating the number of decimal places when rounding the distribution's support.
#' It is necessary to compute posterior probabilities of mode locations.
#' @param tol_mixp Components with a mixture proportion below `tol_mixp` are discarded when estimating modes;
#' note that this does not apply to the biggest component so that it is not possible to discard all components;
#' should be between `0` and `1`; default is `0`.
#' @param tol_x (for continuous mixtures) Tolerance parameter for distance in-between modes; default is `sd(data)/10`
#' where data is the vector of observations from `BayesMix`.
#' If two modes are closer than `tol_x`, only the first estimated mode is kept.
#' @param tol_conv (for continuous mixtures) Tolerance parameter for convergence of the algorithm; default is `1e-8`.
#' @param inside_range Should modes outside of the observations range be discarded? Default is `TRUE`.
#' This sometimes occurs with very small components when K is large.  
#' @return A list of class `bayes_mode` containing:
#'  \item{data}{From `BayesMix`.}
#'  \item{dist}{From `BayesMix`.}
#'  \item{dist_type}{From `BayesMix`.}
#'  \item{pars_names}{From `BayesMix`.}
#'  \item{modes}{Matrix with a row for each draw and columns showing modes.}
#'  \item{p1}{Posterior probability of unimodality.}
#'  \item{tb_nb_modes}{Matrix showing posterior probabilities for the number of modes.}
#'  \item{table_location}{Matrix showing posterior probabilities for mode locations.}
#'  \item{algo}{Algorithm used for mode estimation.}
#'  \item{range}{Range outside which modes are discarded if `inside_range` is `TRUE`.}
#'  \item{BayesMix}{`BayesMix`.}
#' 
#' @details
#' Each draw from the MCMC output after burnin, \eqn{\theta^{(d)}, \quad d = 1,...,D}, leads to a posterior predictive probability
#' density/mass function: 
#' \deqn{p(y | \theta^{(d)}) =\sum_{k=1}^{K} \pi_k^{(d)} p(y | \theta_k^{(d)}).}
#' Using this function, the mode in draw \eqn{d} \eqn{y_{m}^{(d)}}, \eqn{m = 1,..., M^{(d)}},
#' where \eqn{M^{(d)}} is the number of modes, are estimated using the algorithm mentioned
#' in the description above.
#' 
#' After running this procedure across all retained posterior draws, 
#' we compute the posterior probability for the number of modes being \eqn{M} as:
#' \deqn{P(\#\text{modes}=M)=\frac{1}{D}\sum_{d=1}^{D}1(M^{(d)} = M).}
#' Similarly, posterior probabilities for locations of the modes are given by:
#' \deqn{P(y=\text{mode})=\frac{1}{D}\sum_{d=1}^{D} \sum_{m=1}^{M^{(d)}} 1(y = y_m^{(d)}),}
#' for each location \eqn{y} in the range \eqn{[\min(y),\max(y)]}. Obviously,
#' continuous data are not defined on a discrete support;
#' it is therefore necessary to choose a rounding decimal to discretize their support (with the \code{rd} argument).
#' 
#' @references
#' \insertAllCited{}
#
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' 
#' @examples
#' # Example with galaxy data ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = galaxy
#'
#' # estimation
#' bayesmix = bayes_fit(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#' 
#' # mode estimation
#' BayesMode = bayes_mode(bayesmix)
#'
#' # plot 
#' # plot(BayesMode, max_size = 200)
#'
#' # summary 
#' # summary(BayesMode)
#' 
#' # Example with DNA data ================================================
#' set.seed(123) 
#' 
#' # retrieve DNA data
#' y = d4z4
#'
#' # estimation
#' bayesmix = bayes_fit(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "shifted_poisson",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#' 
#' # mode estimation
#' BayesMode = bayes_mode(bayesmix)
#'
#' # plot 
#' # plot(BayesMode, max_size = 200)
#'
#' # summary 
#' # summary(BayesMode)
#' 
#' # Example with a Student t ================================================
#' mu = c(0.5,6)
#' sigma = c(1,2)
#' nu = c(5,5)
#' p = c(0.8,0.2)#'
#' data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
#'          sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
#'
#' fit = c(eta = p, mu = mu, sigma = sigma, nu = nu, xi = c(0,0))
#' fit = rbind(fit, fit)
#' 
#' pdf_func = function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' dist_type = "continuous"
#' 
#' bayesmix = bayes_mixture(fit, data, burnin = 1, 
#' pdf_func = pdf_func, dist_type = dist_type, loc = "mu")
#' 
#' BayesMode = bayes_mode(bayesmix)
#' 
#' # plot 
#' # plot(BayesMode, max_size = 200)
#'
#' # summary 
#' # summary(BayesMode)
#'
#' @export
bayes_mode <- function(BayesMix, rd = 1, tol_mixp = 0, tol_x = sd(BayesMix$data)/10, tol_conv = 1e-8, inside_range = TRUE) {
  assert_that(inherits(BayesMix, "bayes_mixture"), msg = "BayesMix should be an object of class bayes_mixture")
  assert_that(all(c("data", "mcmc", "mcmc_all",
                    "loglik", "K", "dist",
                    "dist_type", "pdf_func", "pars_names",
                    "loc", "nb_var") %in% names(BayesMix)),
              msg = "BayesMix is not a proper bayes_mixture object.") 
  
  assert_that(is.scalar(rd), rd >= 0, round(rd) == rd, msg = "rd should be an integer greater or equal than zero")
  assert_that(is.scalar(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.scalar(tol_mixp) & tol_mixp >= 0 & tol_mixp < 1, msg = "tol_mixp should be a positive scalar between 0 and 1")
  assert_that(is.scalar(tol_conv) & tol_conv > 0, msg = "tol_conv should be a positive scalar")
  
  dist = BayesMix$dist
  data = BayesMix$data
  mcmc = BayesMix$mcmc
  dist_type = BayesMix$dist_type
  pdf_func = BayesMix$pdf_func
  pars_names = BayesMix$pars_names
  loc = BayesMix$loc
  
  if (dist_type == "continuous") {
    range = c(min(data) - sd(data), max(data) + sd(data)) 
  } else {
    range = c(min(data), max(data))
  }
  
  modes = t(apply(mcmc, 1, mix_mode_estimates, dist = dist,
                  pdf_func = pdf_func, dist_type = dist_type,
                  tol_mixp = tol_mixp, tol_x = tol_x, tol_conv = tol_conv,
                  loc = loc, range = range,
                  inside_range = TRUE))

  # Number of modes 
  n_modes = apply(!is.na(modes),1,sum) # number of modes in each MCMC draw
  modes = matrix(modes[, 1:max(n_modes)], nrow = nrow(mcmc))
  colnames(modes) = paste('mode',1:max(n_modes))

  vec_modes = as.vector(modes)
  vec_modes = vec_modes[!is.na(vec_modes)]
  
  if (dist_type == "continuous") {
    if (!is.na(dist) & dist == "normal") {
      algo = "fixed-point"
    } else {
      algo = "Modal Expectation-Maximization (MEM)"
    }
    
    vec_modes = round(vec_modes, rd)
    mode_range = seq(min(vec_modes), max(vec_modes), by = 10^-rd)
  }
  
  if (dist_type == "discrete") {
    mode_range = min(vec_modes):max(vec_modes)
    
    # unique modes to calculate post probs of number of modes
    modes <-  t(apply(mcmc,1,FUN = mix_mode_estimates,
                      range = range,
                      dist = dist,
                      dist_type = dist_type,
                      tol_mixp = tol_mixp,
                      tol_x = tol_x,
                      tol_conv = tol_conv,
                      type = "unique",
                      pdf_func = pdf_func,
                      inside_range = inside_range))
    
    n_modes = apply(!is.na(modes),1,sum)
    
    algo = "discrete"
  }
  
  ### Posterior probability of being a mode for each location
  sum_modes = unlist(lapply(mode_range,
                            FUN = counting,
                            vec = vec_modes))
  
  probs_modes = sum_modes/nrow(mcmc)
  
  table_location = rbind(mode_range, probs_modes)
  rownames(table_location) = c("mode location", "posterior probability")
  
  ##### testing unimodality
  p1 = 0 # posterior probability of unimodality
  
  if(any(n_modes==1)){
    p1 = length(n_modes[n_modes==1])/nrow(modes)
  }
  
  # Test for number of modes : number of modes and their posterior probability
  unique_modes = unique(n_modes) #possible number of modes
  prob_nb_modes = rep(NA_real_,length(unique_modes))
  for (i in 1:length(unique_modes)){
    prob_nb_modes[i] = length(n_modes[n_modes==unique_modes[i]])/nrow(modes)
  }
  tb_nb_modes = rbind(unique_modes,prob_nb_modes)
  rownames(tb_nb_modes) = c("number of modes", "posterior probability")
  
  bayes_mode = list()
  bayes_mode$data = data
  bayes_mode$dist = dist
  bayes_mode$dist_type = dist_type
  bayes_mode$pars_names = pars_names
  bayes_mode$modes = modes
  bayes_mode$p1 = p1
  bayes_mode$tb_nb_modes = tb_nb_modes
  bayes_mode$table_location = table_location
  bayes_mode$algo = algo
  bayes_mode$BayesMix = BayesMix
  bayes_mode$range = range
  
  class(bayes_mode) <- "bayes_mode"
  
  return(bayes_mode)
}

#' @keywords internal
mix_mode_estimates <- function(mcmc, dist = NA_character_, dist_type = NA_character_,
                               tol_mixp, tol_x, tol_conv,
                               pdf_func = NULL, type = "all", range = NULL,
                               loc = NA_character_, inside_range = TRUE) {
  output = rep(NA_real_, length(mcmc))
  
  mix = mixture(mcmc, dist = dist, pdf_func = pdf_func,
                    dist_type = dist_type, range = range, loc = loc)
  modes = mix_mode(mix, tol_mixp, tol_x, tol_conv, type = type)$mode_estimates
  output[1:length(modes)] = modes
  
  return(output)
}

#' @keywords internal
counting <- function(x, vec) {
  length(vec[vec==x])
}