## functions used in the MCMC algorithm

# Posterior of kappa
post_kap <- function(x,LAMBDA,KAPPA) {
  n = length(x) # Number of elements in the component
  result <-  exp(-LAMBDA)*(LAMBDA^(sum(x)-n*KAPPA))/prod(factorial(x-KAPPA)) 
}

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

# Probability mass function for shifted Poisson distribution
spoisspdf <- function(x,LAMDA,KAPPA){
  n = length(x) # number of points
  # Undefined pdf = 0    
  ind1 = which(x<KAPPA)  # spoisspdf not defined for these values
  nundef = length(ind1) # number of undefined points
  pdf1 = matrix(data=0,nrow=nundef,ncol=1) # set pdf = 0 for undefined values
  # Defined pdf        
  ind2 = which(x>=KAPPA) 
  ndef = sum(ind2) # number of undefined points
  x = x[ind2]
  pdf2 = exp(-LAMDA)*(LAMDA^(x-KAPPA))/factorial(x-KAPPA)
  # Combine 
  pdfc = matrix(data=NA,nrow=n,ncol=1)
  pdfc[ind1] = pdf1
  pdfc[ind2] = pdf2
  return(pdfc)           # Return output
}

# Posterior of e0 - Unnormalized target pdf
post_e0 <- function(e0,nu0_p,S0_p,p) {
  K = length(p) # Number of components
  result <-  dgamma(e0,shape = nu0_p, scale = S0_p)*gamma(K*e0)/(gamma(e0)^K)*((prod(p))^(e0-1))
}

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