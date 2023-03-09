#' Bayesian estimation of mixture of skew normal
#' 
#' Gibbs sampler for estimating a SFM of skew normal distributions from xxx 20..
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
#' @importFrom sn dsn
#' @importFrom Rdpack reprompt
#' @importFrom stats median kmeans rgamma rmultinom rnorm 

#' @keywords internal
gibbs_SFM_skew_n <- function(y,
                             K,
                             nb_iter,
                             priors = list(),
                             printing = TRUE){

  # unpacking priors
  # a0 = ifelse(is.null(priors$a0), 10, priors$a0)
  # A0 = ifelse(is.null(priors$A0), a0*K, priors$A0)
  a0 = ifelse(is.null(priors$a0), 1, priors$a0)
  A0 = ifelse(is.null(priors$A0), 200, priors$A0)
  b0 = ifelse(is.null(priors$b0), median(y), priors$b0)
  B0 = ifelse(is.null(priors$B0), (max(y) - min(y))^2, priors$B0)
  c0 = ifelse(is.null(priors$c0), 2.5, priors$c0)
  C0 = ifelse(is.null(priors$C0), 0.5*var(y), priors$C0)
  e0 = ifelse(is.null(priors$e0), a0/A0, priors$e0)
  g0 = ifelse(is.null(priors$g0), 0.5, priors$g0)
  G0 = ifelse(is.null(priors$G0), g0/(0.5*var(y)), priors$G0)
  D_xi = ifelse(is.null(priors$D_xi), 1, priors$D_xi)
  D_psi = ifelse(is.null(priors$D_psi), 1, priors$D_psi)
  
  assert_that(is.scalar(A0) & A0 > 0, msg = "A0 should be positive")

  n_obs = length(y)
  
  #empty objects to store parameters
  sigma2 = rep(NA, K)
  psi = rep(NA, K)
  xi = matrix(NA, nb_iter, K)
  omega = matrix(NA, nb_iter, K)
  alpha = matrix(NA, nb_iter, K)
  eta = matrix(NA, nb_iter, K)
  lp = matrix(0, nb_iter, 1)
  
  # initialisation
  cl_y = kmeans(y, centers = K, nstart = 30)
  
  S <- matrix(0, n_obs, K)
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  b0 = matrix(c(b0, 0),nrow=2)
  B0_inv = diag(1/c(D_xi,D_psi))

  for (k in 1:K){
    xi[1, k] = mean(y[S[, k]==1])
  }
  sigma2 = rep(max(y)^2,K)
  beta = cbind(xi[1, ],0)

  z = matrix(0,n_obs,K)
  cnt_update_e0 = 0
  
  for (m in 1:nb_iter){
    # set base priors
    ck = c0
    Ck = C0
    bk = b0
    Bk = solve(B0_inv)
    
    ## a.1 sample component proportion
    N = colSums(S)
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    probs = matrix(NA, n_obs, K)
    
    for (k in 1:K){
      
      if (N[k] > 0) {
        empty = FALSE
      } else {
        empty = TRUE
      }
      
      #
      
      if(!empty){
        # allocate y
        yk = y[S[ ,k] == 1]
        
        # update z
        Ak = 1/(1 + beta[k, 2]^2/sigma2[k])
        ak = Ak*beta[k, 2]/sigma2[k]*(yk-beta[k, 1])
        # z =  tr_normal(n_obs, ak, diag(Ak, n_obs))
        zk <- MCMCglmm::rtnorm(N[k], ak, sqrt(Ak), lower=0) 

        Xk = matrix(c(rep(1, N[k]), zk), nrow=N[k])
      }
   
      ## a.2 sample xi and psi jointly
      if(!empty){
        Bk = solve(crossprod(Xk)/sigma2[k] + B0_inv)
        bk = Bk%*%(crossprod(Xk, yk)/sigma2[k] + B0_inv%*%b0)
      }
      # beta[k, ] = t(chol(sigma2[k]*Bk)) %*% rnorm(2) + t(bk)
      beta[k, ] = mvtnorm::rmvnorm(1, bk, Bk)
   
      ## a.1 sample Sigma
      if(!empty){
        eps = yk - Xk%*%beta[k, ]
        Ck = C0 + 0.5*sum(eps^2)
        ck = c0 + N[k]/2
      }
      sigma2[k] = 1/rgamma(1, ck, Ck)
      
      # storing
      xi[m, k] = beta[k, 1]
      omega[m, k] = sqrt(sigma2[k] + beta[k, 2]^2)
      alpha[m, k] = beta[k, 2]/sqrt(sigma2[k])
      probs[, k] = eta[m, k] * dsn(y, xi[m, k], omega[m, k], alpha[m, k])
    }
    # browser()
    
    # classification
    pnorm = probs/rowSums(probs)
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1, size = 1, prob = x)))
  
    # 3. sample hyperparameters
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2))
    
    ## SFM: MH step for e0
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
  mcmc_result = cbind(eta, xi, omega, alpha, lp)
  colnames(mcmc_result) = 1:ncol(mcmc_result)
  
  for (i in 1:K){
    colnames(mcmc_result)[c(i, K+i, 2*K+i,3*K+i)] = c(paste0("theta", i),
                                                      paste0("xi", i),
                                                      paste0("omega", i),
                                                      paste0("alpha", i))
  }
  colnames(mcmc_result)[ncol(mcmc_result)] = "loglik"
  
  return(mcmc_result)
}

#' @keywords internal
tr_normal <- function(N,mu,sigma2){
  rd_sample = rnorm(N,mu,sqrt(sigma2))
  rd_sample = rd_sample[rd_sample>=0]
  while (length(rd_sample)<N) {
    rd_sample = c(rnorm(N,mu,sqrt(sigma2)),rd_sample)
    rd_sample = rd_sample[rd_sample>=0] 
  }
  
  return(rd_sample[1:N])
}