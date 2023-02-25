#' Estimation of a mixture of shifted Poisson distributions.
#' 
#' Bayesian estimation of a mixture of shifted Poisson distributions using a Sparse Finite Mixture MCMC algorithm.
#' @param y (a vector of integers) Observations used to fit the model.
#' @param K (an integer) Maximum number of mixture components.
#' @param nb_iter (an integer) Number of MCMC iterations.
#' @param prt print intermediate of the MCMC estimation ? default = TRUE.
#' @returns 
#' A list containing:
#' \itemize{
#'   \item mcmc_draws : a (nb_iter x 3xK) matrix. Parameter draws from the posterior distribution at each MCMC iteration.
#'   \item p_draws : a (nb_iter x K) matrix. Posterior draw of mixture weights at each MCMC iteration.
#'   \item kappa_draws : a (nb_iter x K) matrix. Posterior draws of kappa at each MCMC iteration.
#'   \item lambda_draws : a (nb_iter x K) matrix. Posterior draws of lambda at each MCMC iteration.
#'   \item snj : a (nb_iter x K) matrix. Number of observations allocated to each components at each MCMC iteration.
#'   \item se0 : a vector of size nb_iter. Concentration parameter from (symmetric) Dirichlet distribution at each MCMC iteration.
#'   \item dist : a string indicating the distribution used in the mixture.
#'   \item y, K and nb_iter, given as input.
#' }
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr
#' \insertRef{basturk_bayes_2021}{BayesMultiMode}

#' @importFrom gtools rdirichlet
#' @importFrom stats density dgamma dpois rgamma rmultinom rnorm runif

#' @keywords internal
gibbs_SFM_sp <- function(y, K, nb_iter, prt = TRUE){
  # initialization 
  tmp <- runif(K)
  theta.init <- c(tmp[-K] / sum(tmp), rep(min(y),K), rep(1,K))
  
  names(theta.init) = c(paste("prob ",1:(K-1),sep=""),
                        paste("kappa ",1:K,sep=""),
                        paste("lambda ",1:K,sep=""))
  theta0 <- theta.init
  
  # Error Messages  
  if(round(K) != K | K < 1){
    stop("max number of mixture components should be integer >= 1")
  }
  if(!is.vector(y)){
    stop("data 'y' should be a vector")
  }

  # initial parameters
  ind = K
  p0 <- theta0[1:(K-1)];
  kappa0   <- theta0[ind:(ind+K-1)]; ind=ind+K
  lambda0 <- theta0[ind:(ind+K-1)]; 
  p <- c(p0,(1-sum(p0))); lam <- lambda0; kap <- kappa0
  N <- length(y)     
  lambda.rng <- kappa.rng <- c(0,max(y))
  
  if(any(lambda0<=0))
    stop("'lambda0' should be > 0 for all clusters");
  if(any(kappa0>min(y)))
    stop("'kappa0' should be <= min(y) for all clusters");
  
  ## Priors
  # kappa_j - uniform prior where upper bound varies in MCMC
  kaplb = 0 # Lower bound for kappa
  
  # lambda_j - independent gamma prior - symmetric prior across all components
  nu0 = 5
  invS0 = (nu0-1)
  
  # e0 - independent gamma prior: Eq (3) in MW, FS and Grun
  nu0_p = 1 
  S0_p = 1/200
  
  # Initial conditions
  lam = matrix(data=lam,nrow=1,ncol=K) # lambda 
  kap = matrix(data=kap,nrow=1,ncol=K) # kappa 
  p = matrix(data=p,nrow=1,ncol=K) # probabilities
  Si = t(rmultinom(n = N,size=1,prob=p)) # component indicator
  alpha0 = matrix(data=nu0_p*S0_p,nrow=1,ncol=K) # concentration parameter - note mean of gamma distribution = nu0_p*S0_p
  
  # storage matrices
  slam = matrix(data=NA,nrow=nb_iter,ncol=K) # lambda
  skap = matrix(data=NA,nrow=nb_iter,ncol=K) # kappa
  sP = matrix(data=NA,nrow=nb_iter,ncol=K)   # probabilities
  snj = matrix(data=NA,nrow=nb_iter,ncol=K)  # number of draws in each component
  se0 = matrix(data=NA,nrow=nb_iter,ncol=1) # Symmetric Dirichlet hyperparameter
  sJ = matrix(data=NA,nrow=nb_iter,ncol=1) # Number of non-empty components
  
  # Counters for MH Steps
  cnt_update_kap = matrix(data=0,nrow=1,ncol=K)  # counter for MH step update of kappa
  cnt_update_e0 = 0 # counter for MH step update of concentration parameter
  
  for(m in 1:nb_iter){
    ## Sample kappa, lamda and S for each component, j=1,...,k
    probs = matrix(data=NaN,nrow=N,ncol=K) # Storage for probabilities  
    for (j in 1:K){
      # Set-up
      indj = (Si[,j]==1) # Find all S_i = j
      nj = sum(indj)
      if (nj==0) {
        yj = 0
      } else {
        yj = y[indj]
      }
      
      
      # Sample kappa using MH Step
      kapub = min(yj)
      if (length(kapub) == 0) {# Set to zero if component is empty
        kap[,j] = 0; 
      } else if (kapub < kap[,j]) {# Set to upper bound if outside boudary
        kap[,j] = kapub
      } else {
        temp = draw_kap(yj,lam[,j],kap[,j],kaplb,kapub) # Draw kappa from MH step
        kap[,j] = as.numeric(temp[1])
        cnt_update_kap[,j] = cnt_update_kap[,j] + as.numeric(temp[2])
      }
      
      # Sample lambda from Gamma distribution - change if prior is non-Gamma
      nup = (sum(yj)-nj*kap[,j]) + nu0
      Sp = 1/(nj + invS0)
      lam[,j] = rgamma(1, shape = nup, scale = Sp)
      
      # Sample component indicator from multinomial
      probs[,j] = p[,j]*spoisspdf(y,lam[,j],kap[,j])        
    }
    
    pnorm = probs/replicate(K,rowSums(probs))
    for (jj in 1:N){
      Si[jj,] = t(rmultinom(n = 1,size=1,prob=pnorm[jj,])) # component indicator
    }
    
    ## Sample component probabilities from Dirichlet
    nj2 = colSums(Si) # number of obs in each component
    alphastar = alpha0 + nj2 # posterior concentration parameter
    p = rdirichlet(1,alphastar) 
    
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = alpha0[1] # Symmetric parameter so choose 1 as default
    temp = draw_e0(e0,nu0_p,S0_p,p)
    e0 = as.numeric(temp[1])
    ind_update_e0 = as.numeric(temp[2])
    alpha0 = e0*matrix(data=1,nrow=1,ncol=K) # probabilities
    cnt_update_e0 = cnt_update_e0 + ind_update_e0
    
    ## Ordering based on the expected value (mu = kappa+lambda)
    mu = kap+lam
    matrix_temp = t(rbind(mu,p,kap,lam,nj2))
    matrix_temp = matrix_temp[order(matrix_temp[,1]),]

    ## Store results  
    slam[m,] = matrix_temp[,4]    
    skap[m,] = matrix_temp[,3]     
    sP[m,] = matrix_temp[,2] 
    snj[m,] = matrix_temp[,5] 
    se0[m,] = e0  
    sJ[m,] = K - sum(is.na(slam[m,]))
    
    # Set number of data points in empty components to NaN value
    sP[m, snj[m,]==0] = NA 
    skap[m, snj[m,]==0] = NA 
    slam[m, snj[m,]==0] = NA
    
    colnames(sP) = paste0("theta", 1:ncol(sP))
    colnames(skap) = paste0("kappa", 1:ncol(skap))
    colnames(slam) = paste0("lambda", 1:ncol(slam))
    
    colnames(snj) = paste0("snj", 1:ncol(sP))
    colnames(se0) = "se0"

    ## counter
    if(prt){
      if(m %% (round(nb_iter / 10)) == 0){
        naccept_e0 = 100*cnt_update_e0/m
        cat(paste(100 * m / nb_iter, ' % draws finished. Accept. prob of e0 =', round(naccept_e0), 'percent'), fill=TRUE)
      }
    }
  }
  
  # Store output
  mcmc = cbind(sP, skap, slam, snj, se0)
  
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