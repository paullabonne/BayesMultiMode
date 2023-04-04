#' Bayesian estimation of a mixture of Normal distributions.
#' 
#' MCMC estimation using a sparse finite mixture (SFM) algorithm.
#' 
#' @param y Vector of observations.
#' @param K Maximum number of mixture components.
#' @param nb_iter Number of MCMC iterations.
#' @param priors List of priors. Default is :
#' list(a0 = 1, A0 = 200, b0 = median(y), B0 = (max(y) - min(y))^2, c0 = 2.5, g0 = 0.5, G0 = 100*g0/c0/B0)
#' @param printing Print intermediate output of the MCMC estimation ? default = TRUE.
#' 
#' @returns 
#' mcmc_draws : Parameter draws from the posterior distribution at each MCMC iteration. A (nb_iter x 3K + 1) matrix. 
#' 
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr\cr
#' \insertRef{viallefont2002bayesian}{BayesMultiMode}\cr\cr
#' \insertRef{SFS:Mal:2019}{BayesMultiMode}\cr
#' 
#' @importFrom gtools rdirichlet
#' @importFrom Rdpack reprompt
#' @importFrom stats median kmeans rgamma rmultinom rnorm

#' @keywords internal
gibbs_SFM_normal <- function(y,
                             K,
                             nb_iter,
                             priors = list(),
                             printing = TRUE){
  
  # unpacking priors
  a0 = ifelse(is.null(priors$a0), 1, priors$a0)
  A0 = ifelse(is.null(priors$A0), 200, priors$A0)
  b0 = ifelse(is.null(priors$b0), median(y), priors$b0)
  B0 = ifelse(is.null(priors$B0), (max(y) - min(y))^2, priors$B0)
  c0 = ifelse(is.null(priors$c0), 2.5, priors$c0)
  g0 = ifelse(is.null(priors$g0), 0.5, priors$g0)
  G0 = ifelse(is.null(priors$G0), 100*g0/c0/B0, priors$G0)
  
  # checking priors' validity
  assert_that(is.scalar(A0) & A0 > 0, msg = "A0 should be positive")
  assert_that(is.scalar(B0) & B0 > 0, msg = "B0 should be a positive integer")
  
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
  
  C0 = g0
  e0 = a0/A0
  
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
      
      mu[m, k] = rnorm(1, b, sqrt(B))
      
      # 2. classification
      probs[, k] = eta[m, k] * dnorm(y, mu[m, k], sqrt(sigma2[m, k]))
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) #removed the replicate
    
    ## if the initial classification is bad then some data points won't be 
    # allocated to any components and some rows will be 
    # NAs (because if dividing by zero). We correct this by replacing NAs with
    # equal probabilities
    NA_id = which(is.na(pnorm[,1]))
    pnorm[NA_id, ] = 1/ncol(pnorm)
    
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1,size=1,prob=x)))
    
    # 3. sample hyperparameters
    
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2[m, ]))
    
    ## MH step for e0
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
  mcmc = cbind(eta, mu, sqrt(sigma2), lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K){
    colnames(mcmc)[c(i, K+i, 2*K+i)] = c(paste0("eta", i),
                                         paste0("mu", i),
                                         paste0("sigma", i))
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  
  return(mcmc)
}
