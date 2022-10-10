#' Mode inference using post-processed SFM MCMC draws.
#' 
#' Computes the number of modes, their locations and posterior probabilities.
#' @param theta_draws a (M x 3xJb) matrix. Output of `sfm_mcmc_spmix()` giving MCMC parameter draws after burn-in and discarding empty components.
#' @param y (a vector of integers) Observations used to fit the model.
#' @param mixt (a string) giving the mixture distribution. Default is "shifted_poisson".
#' @returns 
#' A list containing:
#' \itemize{
#'   \item Prob_unimod : Posterior probability of unimodality. (1-Prob_unimod) is equal to the posterior probability of multimodality.
#'   \item table_nb_modes : Possible number of modes and posterior probability for each of those.
#'   \item table_locations : Possible locations of modes and posterior probability for each of those.
#'   \item A list of graphs showing : 
#'   \itemize{
#'        \item 1: The posterior probability of multimodality;
#'        \item 2: Possible number of modes and posterior probability for each of those;
#'        \item 3: Possible locations of modes and posterior probability for each of those;
#' }
#' }
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' # Example with simulated data ================================================
#' #set seed for random number generation
#' set.seed(1) 
#' 
#' # Set the parameters for drawing from a two-component shifted Poisson mixture
#' p1 = 0.3
#' p2 = 1-p1
#' kap1 = 3
#' kap2 = 0
#' lam1 = 1
#' lam2 = 0.5
#' length_data = 70
#' 
#' # Generate data
#' y <- c(rpois(length_data*p1, lam1)+kap1, rpois(length_data*p2, lam2)+kap2)
#' 
#' # Set parameters for the SFM MCMC estimation
#' M = 1000 # Number of MCMC iterations 
#' Jmax = 4 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' 
#' # Proportion of draws burned in
#' S = 0.5
#' 
#' # Post processing
#' post_sfm_mcmc = post_sfm_mcmc(sfm_mcmc, S=S)
#' 
#' # Mode inference
#' bayes_mode(post_sfm_mcmc$theta_draws_slim,y)
#' 
#' # Example with DNA data =====================================================
#' \donttest{
#' y = d4z4
#' M = 5000 # Number of MCMC iterations 
#' Jmax = 10 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' 
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' # Proportion of draws burned in
#' S = 0.5
#' 
#' # Mode inference
#' bayes_mode(post_sfm_mcmc$theta_draws_slim,y)
#' }
#' @export
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{basturk_bayes_2021}{BayesMultiMode}
#' 

bayes_mode <- function(theta_draws, y, mixt="shifted_poisson"){
  Pb <- value <- NULL
  
  if(!is.vector(y)) stop("y should be a vector")
  if(mixt!="shifted_poisson") stop("mixture not suported. Mixtures suported are : shifted_poisson")
  
  
  ### COUNT NUMBER OF MODES
  fn.sub.mixpois<-function(theta.draws_i,y, which.r){
    theta.draws_i = theta.draws_i[!is.na(theta.draws_i)]
    Khat = length(theta.draws_i)/3
    
    if (mixt=="shifted_poisson"){
      theta <- cbind(theta.draws_i[1:Khat],
                     theta.draws_i[(Khat+1):(2*Khat)],
                     theta.draws_i[(2*Khat+1):(3*Khat)])
      
      p <- theta[,1]
      kappa <- theta[,2]
      lambda <-  theta[,3]
      
      ### Getting individual component densities
      pdf.J = matrix(nrow=length(y),ncol=Khat) 
      for(j in 1:Khat){
        pdf.J[,j] = dpois((y-kappa[j]),lambda[j]) * p[j]
      }
    }
    
    ### summing up to get the mixture
    py <- rowSums(pdf.J)
    
    ### Finding the modes
    out <- count_modes(y,py, mode.sel = "leftmode")
    r <- num.modes <- out$num.modes
    if(which.r == 2){
      # mode locations
      r = rep(NA,ncol(theta_draws))
      r[1:length(out$y.peaks)] = out$y.peaks
    }
    if(which.r == 3){
      r = rep(NA,ncol(theta_draws))
      #density at modes
      r[1:length(out$py.peaks)] = out$py.peaks
    }
    if(which.r == 4){
      r = rep(0,length(y))
      
      # account for flat modes
      out_right <- count_modes(y,py, mode.sel = "rightmode")
      for(i in 1:length(out$y.peaks)){
        r[y==out$y.peaks[i]] = 1
        
        if(out$y.peaks[i]!=out_right$y.peaks[i]){
          r[y==out_right$y.peaks[i]] = 1
          
          diff = out_right$y.peaks[i]-out$y.peaks[i]
          if(diff>1){
            for(j in 1:(diff-1))
              r[y==out$y.peaks[i]+j] = 1
          }
        }
      }
      #
    }
    return(r)
  }
  
  # Implement function     
  y.pos <- min(y):max(y) # Range
  n.modes <- apply(theta_draws,1,FUN = fn.sub.mixpois, y = y.pos,which.r = 1) # number modes
  modes <- t(apply(theta_draws,1,FUN = fn.sub.mixpois, y = y.pos,which.r = 2)) # location modes
  modes = as.matrix(modes[, 1:max(n.modes)])
  colnames(modes) = paste('mode',1:max(n.modes))
  p.modes <- t(apply(theta_draws,1,FUN = fn.sub.mixpois, y = y.pos,which.r = 3)) # location modes
  p.modes = as.matrix(p.modes[, 1:max(n.modes)])
  colnames(p.modes) = paste('mode',1:max(n.modes))
  
  # Reshape into a vector
  vec_modes = as.vector(modes)
  vec_modes_without_NAs = vec_modes[!is.na(vec_modes)] # Remove NA values
  
  ## Nb of modes per draw
  nb_modes_per_draw = rep(NA,nrow(modes))
  for(i in 1:nrow(modes)){
    nb_modes_per_draw[i] = length(which(!is.na(modes[i,])))
  }
  
  ##### testing unimodality
  if(any(nb_modes_per_draw==1)){
    Post_prob_number_modes_equal_one = length(n.modes[n.modes==1])/nrow(theta_draws)
  } else {
    Post_prob_number_modes_equal_one = 0
  }
  
  ## Test for number of modes : number of modes and their posterior probability
  possible_nb_modes = unique(nb_modes_per_draw)
  post_prob_nb_modes = rep(NA,length(possible_nb_modes))
  for (i in 1:length(possible_nb_modes)){
    post_prob_nb_modes[i] = length(n.modes[n.modes==possible_nb_modes[i]])/nrow(theta_draws)
  }
  table_nb_modes = rbind(possible_nb_modes,post_prob_nb_modes)
  
  # Number of modes 
  n_modes = apply(!is.na(modes),1,sum) # number of modes in each MCMC draw
  
  
  # Posterior probability of being a mode for each location
  modes_incl_flats <- t(apply(theta_draws,1,FUN = fn.sub.mixpois, y = y.pos, which.r = 4)) # modes including flat ones
  sum_modes_incl_flats = apply(modes_incl_flats,2,sum)
  probs_modes = sum_modes_incl_flats/nrow(theta_draws)
  probs_modes = probs_modes[probs_modes>0]
  range = min(y):max(y)
  location_at_modes = range[sum_modes_incl_flats>0]
  
  table_location = rbind(location_at_modes,probs_modes)
  
  df_g0 = tibble(Pb = "Pb",
                 value = (1-Post_prob_number_modes_equal_one))
  
  g0 = ggplot(data=df_g0, aes(x=Pb, y=value)) +
    ggtitle("Nb. modes > 1") +
    theme_gg + 
    ylim(0, 1) +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity")
  
  df_g1 = as_tibble(t(table_location))
  g1 = ggplot(data=df_g1, aes(x=location_at_modes, y=probs_modes)) +
    theme_gg + 
    ggtitle("Location at the mode") +
    ylim(0, 1) +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity",colour="white")
  
  df_g2 = as_tibble(t(table_nb_modes))
  g2= ggplot(data=df_g2, aes(x=possible_nb_modes, y=post_prob_nb_modes)) +
    theme_gg +
    scale_x_continuous(breaks=possible_nb_modes) +
    ggtitle("Number of modes") +
    ylim(0, 1) +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity")
  
  graphs <- ggarrange(g0, g2, g1,
                      ncol = 3, nrow = 1, widths = c(0.7,1, 1))
  
  return(list(Prob_unimod = Post_prob_number_modes_equal_one,
              table_nb_modes = table_nb_modes,
              table_locations = table_location,
              graphs = graphs))
}
