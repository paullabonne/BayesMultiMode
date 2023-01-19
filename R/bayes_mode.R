#' Bayesian estimation of mixture distributions
#'
#' @param BayesMix object of class `BayesMixture`.
#' @param rd Rounding parameter.
#' @param tol_x Tolerance parameter for small components. Default is 1e-3. All components with mixture weights lower than tol_p are dropped.
#' @param tol_p Tolerance parameter for distance in-between modes. Default is sd(data)/10. If two modes are closer than tol_x, only the first estimated mode is kept.
#' 
#' @return An object of class `BayesMode`.
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' 
#' @export

bayes_mode <- function(BayesMix, rd = 1, tol_p = 1e-3, tol_x = sd(BayesMix$data)/10) {
  assert_that(inherits(BayesMix, "BayesMixture"), msg = "BayesMix should be an object of class BayesMixture")
  assert_that(is.scalar(rd) & rd > 0, msg = "rd should be a positive scalar")
  assert_that(is.vector(tol_p) & tol_p > 0, msg = "tol_p should be a positive scalar")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = "tol_x should be a positive scalar")
  
  dist = BayesMix$dist
  data = BayesMix$data
  mcmc = BayesMix$mcmc
  dist_type = BayesMix$dist_type
  pdf_func = BayesMix$pdf_func
  pars_names = BayesMix$pars_names
    
  assert_that(inherits(BayesMix$mcmc, "draws_matrix"),
              msg = "mcmc in BayesMix is not of type draws_matrix")
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(dist_type %in% c("continuous", "discrete"),
              msg = "dist_type should be either continuous or discrete")
  assert_that(dist %in% c("normal", "student", "skew_t", "poisson",
                          "skew_normal", "shifted_poisson") & is.character(dist),
              msg = "Unsupported distribution. 
              dist should be either normal, student, skew_normal, skew_t, shifted_poisson or NA")
  
  if (dist_type == "continuous") {
    if (dist == "normal") {
      # fixed point
      modes = t(apply(mcmc, 1, fixed_point, data = data, tol_x = sd(data)/10))
    } else {
      # MEM algorithm
      modes = t(apply(mcmc, 1, MEM, dist = dist, data = data, pars_names = pars_names, 
                    pdf_func = pdf_func, tol_x = sd(data)/10, show_plot=F))
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
  }
  
  if (dist == "shifted_poisson") {
    Khat = attr(BayesMix$fit, "K")
    y.pos <- min(data):max(data) # Range
    n.modes <- apply(mcmc,1,FUN = fn.sub.mixpois, y = y.pos, which.r = 1, Khat = Khat) # number modes
    modes <- t(apply(mcmc,1,FUN = fn.sub.mixpois, y = y.pos, which.r = 2, Khat = Khat)) # location modes
    modes = as.matrix(modes[, 1:max(n.modes)], nrow = nrow(mcmc))
    colnames(modes) = paste('mode',1:max(n.modes))

    # Posterior probability of being a mode for each location
    modes_incl_flats <- t(apply(mcmc,1,FUN = fn.sub.mixpois, y = y.pos, which.r = 4, Khat=Khat)) # modes including flat ones
    sum_modes_incl_flats = apply(modes_incl_flats,2,sum)
    probs_modes = sum_modes_incl_flats/nrow(mcmc)
    probs_modes = probs_modes[probs_modes>0]
    range = min(data):max(data)
    location_at_modes = range[sum_modes_incl_flats>0]
    
    table_location = rbind(location_at_modes, probs_modes)
  }
  
  # Number of modes 
  n_modes = apply(!is.na(modes),1,sum) # number of modes in each MCMC draw
  
  ##### testing unimodality
  p1 = 0 #Post_prob_number_modes_equal_one
  
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
