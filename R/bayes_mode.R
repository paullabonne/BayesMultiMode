#' Bayesian mode inference.
#' 
#' Modes are estimated for each mcmc draws and these are used to compute posterior probabilities for the number of modes and their locations.
#' The fixed-point algorithm of Carreira-Perpinan (2000) is used for Gaussian mixtures
#' while the Modal EM algorithm of Li et al. (2007) is used for other continuous mixtures.
#'
#' @param BayesMix object of class `BayesMixture`.
#' @param rd Rounding parameter.
#' @param tol_x Tolerance parameter for distance in-between modes. Default is sd(data)/10.
#' If two modes are closer than tol_x, only the first estimated mode is kept.
#' Not needed for mixtures of discrete distributions.
#' @param show_plot Show density with estimate mode as vertical bars
#' @param nb_iter Number of draws on which the mode-finding algorithm is run. Deault is NULL which means the algorithm is run on all draws.
#' @return An object of class `BayesMode`.
#' 
#'\insertRef{carreira-perpinan_mode-finding_2000}{BayesMultiMode}\cr\cr
#' \insertRef{li_nonparametric_2007}{BayesMultiMode}\cr
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' 
#' @export

bayes_mode <- function(BayesMix, rd = 1, tol_x = sd(BayesMix$data)/10, show_plot = F, nb_iter = NULL) {
  
  assert_that(inherits(BayesMix, "BayesMixture"), msg = "BayesMix should be an object of class BayesMixture")
  assert_that(is.scalar(rd) & rd >= 0, msg = "rd should be greater or equal than zero")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  
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
                          "shifted_poisson", "skew_normal") & is.character(dist),
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
      modes = t(apply(mcmc, 1, fixed_point, data = data, tol_x = tol_x, show_plot = show_plot))
    } else {
      # MEM algorithm
      modes = t(apply(mcmc, 1, MEM, dist = dist, data = data, pars_names = pars_names, 
                      pdf_func = pdf_func, tol_x = tol_x, show_plot = show_plot))
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
                     pdf_func = pdf_func, show_plot = show_plot))
    
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
                      pdf_func = pdf_func, show_plot = show_plot))
    
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
