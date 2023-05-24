#' Bayesian estimation of a mixture of shifted Poisson distributions.
#' 
#' MCMC estimation using a sparse finite mixture (SFM) algorithm.
#' 
#' @param y Vector of discrete observations.
#' @param K Maximum number of mixture components.
#' @param nb_iter Number of MCMC iterations.
#' @param priors List of priors. Default is :
#' list(a0 = 1, A0 = 200, l0 = 1.1, L0 = 1.1/median(y))
#' @param printing Print intermediate output of the MCMC estimation ? default = TRUE.
#' 
#' @returns 
#' mcmc_draws : Parameter draws from the posterior distribution at each MCMC iteration. A (nb_iter x 2K + 1) matrix. 
#' 
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr\cr
#' \insertRef{SFS:Mal:2019}{BayesMultiMode}
#' 
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet
#' @importFrom stats density dgamma dpois rgamma rmultinom rnorm runif

#' @keywords internal
gibbs_SFM_poisson <- function(y,
                              K,
                              nb_iter,
                              priors = list(),
                              printing = TRUE){
  
  # unpacking priors
  a0 = ifelse(is.null(priors$a0), 1, priors$a0)
  A0 = ifelse(is.null(priors$A0), 200, priors$A0)
  l0 = ifelse(is.null(priors$l0), 1.1, priors$l0)
  L0 = ifelse(is.null(priors$L0), 1.1/median(y), priors$L0)
  
  assert_that(is.scalar(A0) & A0 > 0, msg = "A0 should be positive")
  assert_that(is.scalar(L0) & L0 > 0, msg = "L0 should be a positive integer")
  
  # Error Messages  
  if(round(K) != K | K < 1){
    stop("number of mixture components should be integer >= 1")
  }
  
  if(!is.vector(y)){
    stop("data 'y' should be a vector")
  }
  
  assert_that(min(y) > -1, msg = "y should not include negative values")
  
  n_obs <- length(y)
  
  # Initial conditions
  cl_y <- kmeans(y, centers = K, nstart = 30)
  
  S <- matrix(0,length(y),K)
  
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  e0 = a0/A0
  
  # storage matrices
  lambda = matrix(data=NA,nrow=nb_iter,ncol=K) # lambda
  eta = matrix(data=NA,nrow=nb_iter,ncol=K)   # probabilities
  probs = matrix(data=NaN,nrow=n_obs,ncol=K) # Storage for probabilities  
  lp = matrix(0, nb_iter, 1)
  
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
      
      # Sample lambda from Gamma distribution 
      lambda[m,k] = rgamma(1, shape = sum(yk) + l0,
                           rate = N[k] + L0)
      
      # 
      probs[,k] = eta[m,k]*dpois(y,lambda[m,k])
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
    
    ## counter
    if(printing){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc = cbind(eta, lambda, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i)] = c(paste0("eta", i),
                                         paste0("lambda", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  # Return output   
  return(mcmc)
}