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
gibbs_SFM_poisson <- function(y,
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
                           scale = 1/(N[k] + L0))
      
      # 
      probs[,k] = eta[m,k]*dpois(y,lambda[m,k])
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) #removed the replicate
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1,size=1,prob=x)))
    
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = draw_e0(e0,a0,A0,eta[m, ])[[1]]
    
    # compute log lik
    lp[m] = sum(probs)
    
    ## counter
    if(prt){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc = cbind(eta, lambda, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i)] = c(paste0("theta", i),
                                         paste0("lambda", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  # Return output   
  return(mcmc)
}