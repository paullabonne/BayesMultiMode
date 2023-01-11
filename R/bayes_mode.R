#' Bayesian estimation of mixture distributions
#'
#' @param BayesMix object of class `BayesMixture`.
#' 
#' @return An object of class `BayesMode`.
#' 
#' @export

bayes_mode <- function(BayesMix, rd = 1) {
  stopifnot(class(BayesMix)=="BayesMixture")
  
  dist = BayesMix$dist
  data = BayesMix$data
  mcmc = BayesMix$mcmc
  dist_type = BayesMix$dist_type
    
  BayesMode = list()
  BayesMode$data = data
  BayesMode$dist = dist
  BayesMode$dist_type = dist_type
  
  if (dist_type == "continuous") {
    if (dist == "normal") {
      ### fixed point
      modes = t(apply(mcmc, 1, fixed_point, y = data, tol_x = sd(data)/10)) 
    } 
    
    if (dist %in% c("student", "skew_normal")) {
      ### MEM 
      modes = t(apply(mcmc, 1, MEM, dist = dist, y = data, tol_x = sd(data)/10, show_plot=F))
    }
    
    ### Posterior probability of being a mode for each location
    m_range = seq(from = min(round(data,rd)), to = max(round(data,rd)), by = 1/(10^rd)) # range of potential values for the modes
    modes_disc = round(modes, digits = rd)
    
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
    modes <- t(apply(mcmc,1,FUN = fn.sub.mixpois, y = y.pos,which.r = 2, Khat = Khat)) # location modes
    modes = as.matrix(modes[, 1:max(n.modes)])
    colnames(modes) = paste('mode',1:max(n.modes))
    BayesMode$modes = modes
    
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
