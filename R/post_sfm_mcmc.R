#' Simple post-processing of SFM MCMC output.
#' 
#' Gives a matrix of MCMC parameters after burn-in and discarding empty components.
#' @param sfm_mcmc a list. Output of `sfm_mcmc_spmix()` containing the parameter draws from the posterior distribution at each MCMC iteration.
#' @param S (a number between 0 and 1) The first S*M draws will be discarded as a burn-in. M is the total number of MCMC iterations.
#' @returns 
#' A list containing:
#' \itemize{
#'   \item A (M x 3xJb) matrix. Returns theta_draws after burn-in (discarding) the S*M rows. M is the number of rows of theta_draws (number of MCMC iterations). Jb is the number of components which are non-empty in at least one of the draws.
#'   \item J_ne (an integer). The number of non-empty components.
#' }
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
#' # Post processing
#' post_sfm_mcmc = post_sfm_mcmc(sfm_mcmc, S=S)
#' }
#' @export

post_sfm_mcmc <- function(sfm_mcmc, S){
  if(S<=0 || S>1){
    stop("S should be number greater than 0 and inferior or equal to 1")
  }
  
  #Number of draws
  M = nrow(sfm_mcmc$theta_draws)

  #Number of draws to be discarded
  Nb_Burned = M*S

  #Number of components including empty ones
  Jmax = sfm_mcmc$Jmax
  
  #Burn in
  theta_draws_slim = sfm_mcmc$theta_draws[(Nb_Burned+1):M,]
  
  # Discard empty components
  # when a component is empty in a given draw it has a NA. 
  # Thus we discard components we are always NAs (empty in all draws).
  temp = theta_draws_slim
  temp[!is.na(temp)] = 1
  temp[is.na(temp)] = 0
  sum_temp = apply(temp,2,sum)
  
  theta_draws_slim = theta_draws_slim[,sum_temp>0]
  
  # number of non-tempty components
  if(sfm_mcmc$mixt=="shifted_poisson"){
    J_ne = length(which(sum_temp>0))/3
  }

  return(list(theta_draws_slim = theta_draws_slim,
         J_ne = J_ne))
}
