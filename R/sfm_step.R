## functions used in the SFM MCMC algorithm

#' @keywords internal
# Posterior of kappa
post_kap <- function(x,LAMBDA,KAPPA) {
  n = length(x) # Number of elements in the component
  result <-  exp(-LAMBDA)*(LAMBDA^(sum(x)-n*KAPPA))/prod(factorial(x-KAPPA)) 
}

#' @keywords internal
# Draw kappa from posterior using MH step
draw_kap <- function(x,LAMBDA,KAPPA,kaplb,kapub) {
  n = length(x) # Number of elements in the component
  KAPPAhat = sample(kaplb:kapub, 1) # Draw candidate from uniform proposal
  accratio = post_kap(x,LAMBDA,KAPPAhat)/post_kap(x,LAMBDA,KAPPA) # Acceptance ratio
  if(is.na(accratio)){
    accprob = 1 # Acceptance probability if denominator = inf (numerical error due to large number) 
  } else {
    accprob = min(c(accratio,1)) # Acceptance probability if denominator not 0
  }
  rand = runif(1, min = 0, max = 1) # Random draw from unifrom (0,1)
  if (rand < accprob) {
    KAPPA = KAPPAhat; # update
    acc = 1; # indicate update
  } else{
    KAPPA = KAPPA; # don't update
    acc = 0;  # indicate failure to update
  }
  out <- list(KAPPA, acc) # Store output in list
  return(out)           # Return output
}

# Posterior of e0 - Unnormalized target pdf
post_e0 <- function(e0,nu0_p,S0_p,p) {
  K = length(p) # Number of components
  result <-  dgamma(e0,shape = nu0_p, scale = S0_p)*gamma(K*e0)/(gamma(e0)^K)*((prod(p))^(e0-1))
}

#' @keywords internal
# Draw from the unnormalized target pdf for hyperparameter e0
draw_e0 <- function(e0,nu0,S0,p){
  # Define terms
  e0hat = e0 + rnorm(1,0,0.1) # Draw a candidate from Random walk proposal
  accratio = post_e0(e0hat,nu0,S0,p)/post_e0(e0,nu0,S0,p) # Acceptance ratio
  if(is.na(accratio)){
    accprob = 0 # Acceptance probability if denominator = 0 
  } else {
    accprob = min(c(accratio,1)) # Acceptance probability if denominator not 0
  }
  # MH Step
  rand = runif(1, min = 0, max = 1) # Random draw from unifrom (0,1)
  if (rand < accprob) {
    e0 = e0hat; # update
    acc = 1; # indicate update
  } else{
    e0 = e0; # don't update
    acc = 0;  # indicate failure to update
  }
  out <- list(e0, acc) # Store output in list
  return(out)           # Return output
}