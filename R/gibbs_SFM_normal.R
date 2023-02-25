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
                            e0 = 0.01,
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
  # S <- rbind(tabulate(cl_y$cluster),
  #            matrix(0,length(y)-1,K))
  mu[1,] <- cbind(t(cl_y$centers))
  eta[1,] = rep(1/K, K)
  S = t(rmultinom(n = length(y),size=1,prob=eta[1,])) # component indicaton
  C0 = g0 #not sure
  alpha0 = rep(a0*A0, K)
  cnt_update_e0 = 0 # counter for MH step update of concentration parameter
  
  # sampling
  for (m in 2:nb_iter){
    # 1. parameter simulation conditional on the classification
    
    ## a. sample component proportion
    N = colSums(S)
    
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    probs = matrix(NA, length(y), K)
    for (k in 1:K){
      
      ## b. sample sigma
      v = y[S[, k]==1] - mu[m-1, k]
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
      # S = apply(probs, 1, function(x) sample(1:K, 1, prob = x))
    }
    
    pnorm = probs#/rowSums(probs) #removed the replicate
    for (jj in 1:length(y)){
      S[jj,] = t(rmultinom(n = 1,size=1,prob=pnorm[jj,])) # component indicator
    }
    
    # 3. sample hyperparameters
    
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2[m, ]))
    ## b. MH step for e0
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    MH_step = draw_e0(e0[1],a0,A0,eta[m, ])
    alpha0 = as.numeric(MH_step[1])
    ind_update_e0 = as.numeric(MH_step[2])
    e0 = alpha0*matrix(data=1,nrow=1,ncol=K) # probabilities
    cnt_update_e0 = cnt_update_e0 + ind_update_e0
    
    ## c. sample lambda, the penalty parameter for b0 / skip for now
    
    ## d. sample b0 / skip for now
    
    # 4. ramdom permutation of the labeling / skip
    
    # compute log lik
    lp[m] = sum(probs)
    
    ## counter
    if(prt){
      if(m %% (round(nb_iter / 10)) == 0){
        naccept_e0 = 100*cnt_update_e0/m
        cat(paste(100 * m / nb_iter, ' % draws finished. Accept. prob of e0 =', round(naccept_e0), 'percent'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc_result = cbind(eta, mu, sqrt(sigma2), lp)
  colnames(mcmc_result) = 1:ncol(mcmc_result)
  
  for (i in 1:K){
    colnames(mcmc_result)[c(i, K+i, 2*K+i)] = c(paste0("theta", i),
                                                paste0("mu", i),
                                                paste0("sigma", i))
  }
  colnames(mcmc_result)[ncol(mcmc_result)] = "loglik"
  
  return(mcmc_result)
}
