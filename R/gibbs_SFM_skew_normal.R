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
                             a0 = 1,
                             A0 = 1/200,
                             b0 = median(y),
                             c0 = 2.5,
                             C0 = 0.5*var(y),
                             e0 = a0*A0,
                             g0 = 0.5,
                             G0 = g0/(0.5*var(y)),
                             D_xipsi = 10,
                             prt = TRUE){
  
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
  
  # b0 = matrix(c(b0, 0),nrow=2)
  b0 = matrix(c(0, 0),nrow=2)
  B0_inv = diag(1/c(D_xipsi,D_xipsi))

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
      
      # b.1.1 classification
      xi[m, k] = beta[k, 1]
      omega[m, k] = sqrt(sigma2[k] + beta[k, 2]^2)
      alpha[m, k] = beta[k, 2]/sqrt(sigma2[k])
      probs[, k] = eta[m, k] * dsn(y, xi[m, k], omega[m, k], alpha[m, k])
    }
    # browser()
    
    # b.1.2 classification
    pnorm = probs/rowSums(probs)
    S = t(apply(pnorm, 1, function(x) rmultinom(n = 1, size = 1, prob = x)))
  
    # 3. sample hyperparameters
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2))
    
    ## SFM: MH step for e0
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    MH_step = draw_e0(e0,a0,A0,eta[m, ])
    e0 = MH_step[[1]]
    ind_update_e0 = MH_step[[2]]
    cnt_update_e0 = cnt_update_e0 + ind_update_e0

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

#' @keywords snbis
SN_bis <- function(y,xi,omega,alpha){
  2*dnorm(y-xi,sd=sqrt(omega))*pnorm(alpha/sqrt(omega)*(y-xi))
}