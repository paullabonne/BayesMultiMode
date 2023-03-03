#' Bayesian estimation of mixture of normals
#' 
#' Gibbs sampler for estimating a SFM of normal distributions from xxx 20.. (univariate version)
#' @param y (a vector of integers) Observations used to fit the model.
#' @param K (an integer) Maximum number of mixture components.
#' @param nb_iter (an integer) Number of MCMC iterations.
#' @param prt print intermediate of the MCMC estimation ? default = TRUE.
#' @returns 
#' A matrix of MCMC samples
#' 
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr
#' 
#' @importFrom gtools rdirichlet
#' @importFrom Rdpack reprompt
#' @importFrom stats median kmeans rgamma rmultinom rnorm

#' @keywords internal
gibbs_SFM_normal <- function(y,
                            K,
                            nb_iter,
                            a0 = 1,
                            A0 = 1/200,
                            b0 = median(y),
                            B0 = (max(y) - min(y))^2,
                            c0 = 2.5,
                            e0 = a0*A0,
                            g0 = 0.5,
                            G0 = 100*g0/c0/B0,
                            prt = TRUE){
   
  #empty objects to store parameters
  mu = matrix(NA, nb_iter, K)
  sigma2 = matrix(NA, nb_iter, K)
  eta = matrix(NA, nb_iter, K)
  lp = matrix(0, nb_iter, 1)
  
  # initialisation
  cl_y <- kmeans(y, centers = K, nstart = 30)
  
  S <- matrix(0,length(y),K)
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  mu[1,] <- cbind(t(cl_y$centers))

  C0 = g0 #not sure

  # sampling
  for (m in 2:nb_iter){
    # 1. parameter simulation conditional on the classification
    
    ## a. sample component proportion
    N = colSums(S)
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    probs = matrix(NA, length(y), K)
    for (k in 1:K){
      
      ## b. sample sigma
      ck = c0 + N[k]/2
      Ck = C0 + 0.5*sum((y[S[, k]==1]-mu[m-1, k])^2)
      sigma2[m, k] = 1/rgamma(1, ck, Ck)
      
      ## c. sample mu
      B = 1/(1/B0 + N[k]/sigma2[m, k])
      if (N[k]!=0){
        b = B*(b0/B0 + N[k]*mean(y[S[, k]==1])/sigma2[m, k])
      } else {
        b = B*(b0/B0)
      }
      
      mu[m, k] = rnorm(1, b, B)
      
      # 2. classification
      probs[, k] = eta[m, k] * dnorm(y, mu[m, k], sqrt(sigma2[m, k]))
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) #removed the replicate
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1,size=1,prob=x)))
    
    # 3. sample hyperparameters
    
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2[m, ]))
    
    ## MH step for e0
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = draw_e0(e0,a0,A0,eta[m, ])[[1]]

    ## c. sample lambda, the penalty parameter for b0 / skip for now
    
    ## d. sample b0 / skip for now
    
    # 4. ramdom permutation of the labeling / skip
    
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
  mcmc = cbind(eta, mu, sqrt(sigma2), lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i, 2*K+i)] = c(paste0("theta", i),
                                                paste0("mu", i),
                                                paste0("sigma", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  return(mcmc)
}
