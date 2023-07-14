#' Bayesian mode inference
#' 
#' This function estimates modes for each mcmc draw and uses these estimates to compute posterior
#' probabilities for the number of modes and their locations (following the approach of Cross et al. 2023).
#' The fixed-point algorithm of Carreira-Perpinan (2000) is used for Gaussian mixtures.
#' The Modal EM algorithm of Li et al. (2007) is used for other continuous mixtures.
#' A basic algorithm is used for discrete mixtures (see Cross et al. 2023).
#'
#' @param BayesMix An object of class \code{BayesMixture}.
#' @param rd Rounding parameter.
#' @param tol_x Tolerance parameter for distance in-between modes; default is sd(data)/10 where data is an element of argument \code{BayesMix}.
#' If two modes are closer than \code{tol_x}, only the first estimated mode is kept.
#' @param tol_conv Tolerance parameter for convergence of the algorithm; default is 1e-8.
#' Not needed for mixtures of discrete distributions.
#' @param show_plot Show density with estimated mode as vertical bars ?
#' @param nb_iter Number of draws on which the mode-finding algorithm is run; default is NULL which means the algorithm is run on all draws.
#' @return A list of class \code{BayesMode} containing
#' \itemize{
#'  \item{data}{ - from \code{BayesMix}.}
#'  \item{dist}{ - from \code{BayesMix}.}
#'  \item{dist_type}{ - from \code{BayesMix}.}
#'  \item{pars_names}{ - from \code{BayesMix}.}
#'  \item{modes}{ - Matrix with a row for each draw and columns showing modes.}
#'  \item{p1}{ - Posterior probability of unimodality.}
#'  \item{tb_nb_modes}{ - Matrix showing posterior probabilities for the number of modes.}
#'  \item{table_location}{ - Matrix showing the posterior probabilities for location points being modes.}
#' }
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
#' \insertRef{cross_2023}{BayesMultiMode}
#' \insertRef{carreira-perpinan_mode-finding_2000}{BayesMultiMode}\cr\cr
#' \insertRef{li_nonparametric_2007}{BayesMultiMode}
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
#' bayesmix = bayes_estimation(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#' 
#' # mode estimation
#' bayesmode = bayes_mode(bayesmix)
#'
#' # plot 
#' # plot(bayesmode, max_size = 200)
#'
#' # summary 
#' # summary(bayesmode)
#' 
#' # Example with DNA data ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = d4z4
#'
#' # estimation
#' bayesmix = bayes_estimation(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "shifted_poisson",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#' 
#' # mode estimation
#' bayesmode = bayes_mode(bayesmix)
#'
#' # plot 
#' # plot(bayesmode, max_size = 200)
#'
#' # summary 
#' # summary(bayesmode)
#' 
#' # Example with a Student t ================================================
#' mu = c(0.5,6)
#' sigma = c(1,2)
#' nu = c(5,5)
#' p = c(0.8,0.2)
#' params = c(eta = p, mu = mu, sigma = sigma, nu = nu)
#' pars_names = c("eta", "mu", "sigma", "nu")
#' dist_type = "continuous"
#'
#' data = c(sn::rst(p[1]*1000, mu[1], sigma[1], nu = nu[1]),
#'          sn::rst(p[2]*1000, mu[2], sigma[2], nu = nu[2]))
#'
#' fit = c(eta = p, mu = mu, sigma = sigma, nu = nu)
#' fit = rbind(fit, fit)
#' 
#' pdf_func = function(x, pars) {
#'   sn::dst(x, pars["mu"], pars["sigma"], pars["xi"], pars["nu"])
#' }
#' 
#' bayesmix = new_BayesMixture(fit, data, K = 2, burnin = 1,
#' pars_names = pars_names, pdf_func = pdf_func, dist_type = dist_type)
#' 
#' bayesmode = bayes_mode(bayesmix)
#' 
#' # plot 
#' # plot(bayesmode, max_size = 200)
#'
#' # summary 
#' # summary(bayesmode)
#'
#' @export

bayes_mode <- function(BayesMix, rd = 1, tol_x = sd(BayesMix$data)/10, tol_conv = 1e-8, show_plot = FALSE, nb_iter = NULL) {
  assert_that(inherits(BayesMix, "BayesMixture"), msg = "BayesMix should be an object of class BayesMixture")
  assert_that(is.scalar(rd) & rd >= 0, msg = "rd should be greater or equal than zero")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  assert_that(is.logical(show_plot), msg = "show_plot should be either TRUE or FALSE")
  
  if (!is.null(nb_iter)) {
    assert_that(is.scalar(nb_iter) & nb_iter > 0, msg = "nb_iter should be a positive integer") 
  }

  dist = BayesMix$dist
  data = BayesMix$data
  mcmc = BayesMix$mcmc
  dist_type = BayesMix$dist_type
  pdf_func = BayesMix$pdf_func
  pars_names = BayesMix$pars_names
  
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(dist_type %in% c("continuous", "discrete"),
              msg = "dist_type should be either continuous or discrete")
  assert_that(dist %in% c("normal", "poisson",
                          "shifted_poisson", "skew_normal", "NA") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, skew_normal, poisson,
              shifted_poisson, or NA")
  
  
  # if nb_iter is specified (find the mode on a limited number of iterations):
  if (!is.null(nb_iter)) {
    mcmc = mcmc[sample(1:nrow(mcmc), nb_iter), , drop = FALSE]
  }
  
  if (dist_type == "continuous") {
    if (dist == "normal") {
      # fixed point
      modes = t(apply(mcmc, 1, fixed_point, data = data, pars_names = pars_names,
                      tol_x = tol_x, tol_conv = tol_conv, show_plot = show_plot))
    } else {
      # MEM algorithm
      modes = t(apply(mcmc, 1, MEM, dist = dist, data = data, pars_names = pars_names, 
                      pdf_func = pdf_func, tol_x = tol_x, tol_conv = tol_conv, show_plot = show_plot))
    }
    
    ### Posterior probability of being a mode for each location
    m_range = seq(from = min(round(data,rd)), to = max(round(data,rd)), by = 1/(10^rd)) # range of potential values for the modes
    modes_disc = round(modes, rd)
    
    matrix_modes = matrix(0, nrow = nrow(modes), ncol = length(m_range))
    for (i in 1:nrow(matrix_modes)) {
      matrix_modes[i, modes_disc[i, ][!is.na(modes_disc[i, ])] %.in% m_range] = 1
    }
    
    sum_modes = apply(matrix_modes,2,sum)
    probs_modes = sum_modes/nrow(modes)
    probs_modes = probs_modes[probs_modes>0]
    location_at_modes = m_range[sum_modes>0]
    table_location = rbind(location_at_modes, probs_modes)
    
    # Number of modes 
    n_modes = apply(!is.na(modes),1,sum) # number of modes in each MCMC draw
  }
  
  if (dist_type == "discrete") {

    # Posterior probability of being a mode for each location
    modes <- t(apply(mcmc,1,FUN = discrete_MF, data = data,
                     pars_names = pars_names, dist = dist,
                     pmf_func = pdf_func, show_plot = show_plot))
    
    modes_xid = matrix(0, nrow(modes), ncol(modes))
    x = min(data):max(data)
    for (i in 1:nrow(modes)) {
      modes_xid[i, which(x %in% na.omit(modes[i,]))] = 1
    }
 
    # number of modes
    n_modes = rowSums(modes_xid)
    
    sum_modes = apply(modes_xid,2,sum)
    probs_modes = sum_modes/nrow(mcmc)
    probs_modes = probs_modes[probs_modes>0]
    x = min(data):max(data)
    location_at_modes = x[sum_modes>0]
    
    modes = as.matrix(modes[, 1:max(n_modes)], nrow = nrow(mcmc))
    colnames(modes) = paste('mode',1:max(n_modes))
    
    table_location = rbind(location_at_modes, probs_modes)
    
    # unique modes to calculate post probs of number of modes
    modes <-  t(apply(mcmc,1,FUN = discrete_MF, data = data, type = "unique",
                      pars_names = pars_names, dist = dist,
                      pmf_func = pdf_func, show_plot = show_plot))
    
    n_modes = apply(!is.na(modes),1,sum)
  }
  
  ##### testing unimodality
  p1 = 0 # posterior probability of unimodality
  
  if(any(n_modes==1)){
    p1 = length(n_modes[n_modes==1])/nrow(modes)
  }
  
  # Test for number of modes : number of modes and their posterior probability
  unique_modes = unique(n_modes) #possible number of modes
  prob_nb_modes = rep(NA,length(unique_modes))
  for (i in 1:length(unique_modes)){
    prob_nb_modes[i] = length(n_modes[n_modes==unique_modes[i]])/nrow(modes)
  }
  tb_nb_modes = rbind(unique_modes,prob_nb_modes)
  
  BayesMode = list()
  BayesMode$data = data
  BayesMode$dist = dist
  BayesMode$dist_type = dist_type
  BayesMode$pars_names = pars_names
  BayesMode$modes = modes
  BayesMode$p1 = p1
  BayesMode$tb_nb_modes = tb_nb_modes
  BayesMode$table_location = table_location
  
  class(BayesMode) <- "BayesMode"
  
  return(BayesMode)
}

## This function overcomes the problem arising from comparing floating points
## (0.1==0.1 can be false for instance,
## see https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal/9508558#9508558)
#' @keywords internal
`%.in%` = function(a, b, eps = sqrt(.Machine$double.eps)) {
  output = rep(F,length(b))
  for (x in a){
    output <- (abs(b-x) <= eps) | output
  }
  output
}
