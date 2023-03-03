#' Estimation of a mixture of Poisson distributions.
#' 
#' Bayesian estimation of a mixture of Poisson distributions using a Sparse Finite Mixture MCMC algorithm.
#' @param y (a vector of integers) Observations used to fit the model.
#' @param K (an integer) Maximum number of mixture components.
#' @param nb_iter (an integer) Number of MCMC iterations.
#' @param prt print intermediate of the MCMC estimation ? default = TRUE.
#' 
#' @returns 
#' mcmc_draws : a (nb_iter x 2K + 1) matrix. Parameter draws from the posterior distribution at each MCMC iteration.
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr
#' \insertRef{basturk_bayes_2021}{BayesMultiMode}

#' @importFrom gtools rdirichlet
#' @importFrom stats density dgamma dpois rgamma rmultinom rnorm runif

#' @keywords internal
gibbs_SFM_sp <- function(y,
                         K,
                         nb_iter,
                         a0 = 1,
                         A0 = 1/200,
                         l0 = 5,
                         L0 = l0 - 1,
                         e0 = a0*A0,
                         prt = TRUE){
  
  # Error Messages  
  if(round(K) != K | K < 1){
    stop("number of mixture components should be integer >= 1")
  }
  
  if(!is.vector(y)){
    stop("data 'y' should be a vector")
  }
  
  n_obs <- length(y)
  
  # Initial conditions
  cl_y <- kmeans(y, centers = K, nstart = 30)
  
  S <- matrix(0,length(y),K)
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  # storage matrices
  kappa = matrix(data=NA,nrow=nb_iter,ncol=K) 
  lambda = matrix(data=NA,nrow=nb_iter,ncol=K) 
  eta = matrix(data=NA,nrow=nb_iter,ncol=K)   
  probs = matrix(data=NaN,nrow=n_obs,ncol=K)
  lp = matrix(0, nb_iter, 1)
  
  kappa_m = rep(0,K)
  lambda_m = rep(1,K)
  
  ## Sample lamda and S for each component, k=1,...,k
  for(m in 1:nb_iter){
    
    # Compute number of observations allocated in each component
    N = colSums(S)
    
    ## sample component proportion
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    for (k in 1:K){
      
      if (N[k]==0) {
        yk = 0
      } else {
        yk = y[S[, k]==1]
      }
      
      # Sample kappa using MH Step
      kapub = min(yk)
      if (length(kapub) == 0) {# Set to zero if component is empty
        kappa_m[k] = 0; 
      } else if (kapub < kappa_m[k]) {# Set to upper bound if outside boudary
        kappa_m[k] = kapub
      } else {
        temp = draw_kap(yk, lambda_m[k], kappa_m[k], kaplb = 0, kapub) # Draw kappa from MH step
        kappa_m[k] = temp[[1]]
      }
      
      # Sample lambda from Gamma distribution 
      lambda_m[k] = rgamma(1, shape = sum(yk) - N[k]*kappa_m[k] + l0,
                           scale = 1/(N[k] + L0))
      
      # 
      probs[,k] = eta[m,k]*dpois(y - kappa_m[k], lambda_m[k])
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) #removed the replicate
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1,size=1,prob=x)))
    
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = draw_e0(e0,a0,A0,eta[m, ])[[1]]
    
    # compute log lik
    lp[m] = sum(probs)
    
    # storing
    lambda[m,] = lambda_m
    kappa[m, ] = kappa_m
    
    ## counter
    if(prt){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc = cbind(eta, kappa, lambda, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i, 2*K+i)] = c(paste0("theta", i),
                                         paste0("kappa", i),
                                         paste0("lambda", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  # Return output   
  return(mcmc)
}

## functions used in the MCMC algorithm

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

#' @keywords internal
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