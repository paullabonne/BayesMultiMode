#' Bayesian estimation of a mixture of Poisson distributions.
#' 
#' MCMC estimation using a sparse finite mixture (SFM) algorithm.
#' 
#' @param y Vector of discrete observations.
#' @param K Maximum number of mixture components.
#' @param nb_iter Number of MCMC iterations.
#' @param priors List of priors. Default is :
#' list(a0 = 1, A0 = 200, l0 = 5, L0 = l0 - 1)
#' @param printing Print intermediate output of the MCMC estimation ? default = TRUE.
#' 
#' @returns 
#' mcmc_draws : Parameter draws from the posterior distribution at each MCMC iteration. A (nb_iter x 3K + 1) matrix. 
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr\cr
#' \insertRef{SFS:Mal:2019}{BayesMultiMode}
#' 
#' @importFrom gtools rdirichlet
#' @importFrom stats density dgamma dpois rgamma rmultinom rnorm runif

#' @keywords internal
gibbs_SFM_sp <- function(y,
                         K,
                         nb_iter,
                         priors = list(),
                         printing = TRUE){
  
  # unpacking priors
  a0 = ifelse(is.null(priors$a0), 1, priors$a0)
  A0 = ifelse(is.null(priors$A0), 200, priors$A0)
  l0 = ifelse(is.null(priors$l0), 5, priors$l0)
  L0 = ifelse(is.null(priors$L0), l0 - 1, priors$L0)
  
  assert_that(is.scalar(A0) & A0 > 0, msg = "A0 should be positive")
  assert_that(is.scalar(L0) & L0 > 0, msg = "L0 should be a positive integer")
  
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
  
  e0 = a0/A0
  
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
    pnorm = probs/rowSums(probs)
    
    ## if the initial classification is bad then some data points won't be 
    # allocated to any components and some rows will be 
    # NAs (because if dividing by zero). We correct this by replacing NAs with
    # equal probabilities
    NA_id = which(is.na(pnorm[,1]))
    pnorm[NA_id, ] = 1/ncol(pnorm)
    
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1,size=1,prob=x)))
    
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = draw_e0(e0,a0,1/A0,eta[m, ])[[1]]
    
    # compute log lik
    lp[m] = sum(probs)
    
    # storing
    lambda[m,] = lambda_m
    kappa[m, ] = kappa_m
    
    ## counter
    if(printing){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc = cbind(eta, kappa, lambda, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i, 2*K+i)] = c(paste0("eta", i),
                                         paste0("kappa", i),
                                         paste0("lambda", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  # Return output   
  return(mcmc)
}