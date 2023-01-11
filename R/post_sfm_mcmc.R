#' Simple post-processing of SFM MCMC output.
#' 
#' Gives a matrix of MCMC parameters after burn-in and discarding empty components.
#' @param sfm_mcmc a list. Output of `sfm_mcmc_spmix()` containing the parameter draws from the posterior distribution at each MCMC iteration.
#' @param burnin (positive integer) Number of draws used for burnin
#' @returns 
#' A list containing:
#' \itemize{
#'   \item A (M x 3xJb) matrix. Returns theta_draws after burn-in (discarding) the S*M rows. M is the number of rows of theta_draws (number of MCMC iterations). Jb is the number of components which are non-empty in at least one of the draws.
#' }
#' 
#' @export

post_sfm_mcmc <- function(mcmc){
  
  mcmc_post = mcmc
  
  #Number of draws
  M = nrow(mcmc)
  
  mcmc_par <- attributes(mcmc)
  # browser()
  #Burn in
  mcmc_post = mcmc_post[(mcmc_par$warmup+1):M,]
  
  # Discard empty components
  # when a component is empty in a given draw it has a NA. 
  # Thus we discard components we are always NAs (empty in all draws).
  temp = mcmc_post
  temp[!is.na(temp)] = 1
  temp[is.na(temp)] = 0
  sum_temp = apply(temp,2,sum)
  
  mcmc_post = mcmc_post[,sum_temp>0]
  
  # number of non-tempty components
  # if(sfm_mcmc$mixt=="shifted_poisson"){
  #   K_non_emp = length(which(sum_temp>0))/3
  # }

  return(mcmc_post)
}
