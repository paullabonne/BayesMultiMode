#' SFM MCMC algorithms to estimate mixtures.
#' 
#' @importFrom gtools rdirichlet
#' @importFrom Rdpack reprompt
#' @importFrom stats median kmeans rgamma rmultinom rnorm density dgamma dpois runif
#' @importFrom MCMCglmm rtnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom sn dsn

# wrapper
#' @keywords internal
gibbs_SFM <- function(y,
                      K,
                      nb_iter,
                      nb_burnin,
                      priors = list(),
                      print = TRUE,
                      dist) {
  
  if (dist == "normal") {
    mcmc = gibbs_SFM_normal(y, K, nb_iter, priors, print)
  }
  
  if (dist == "skew_normal") {
    mcmc = gibbs_SFM_skew_n(y, K, nb_iter, priors, print)
  }
  
  if (dist == "poisson") {
    mcmc = gibbs_SFM_poisson(y, K, nb_iter, priors, print)
  }
  
  if (dist == "shifted_poisson") {
    mcmc = gibbs_SFM_shifted_poisson(y, K, nb_iter, priors, print)
  }

  if (dist == "neg_binomial") {
    mcmc = gibbs_SFM_neg_binomial(y, K, nb_iter, nb_burnin, priors, print)
  }

  if (dist == "shifted_neg_binomial") {
    mcmc = gibbs_SFM_shifted_neg_binomial(y, K, nb_iter, nb_burnin, priors, print)
  }
  
  return(mcmc)
  
}


#' @keywords internal
gibbs_SFM_normal <- function(y,
                             K,
                             nb_iter,
                             priors = list(),
                             print = TRUE){
  
  # unpacking priors
  a0 = priors$a0
  A0 = priors$A0
  b0 = priors$b0
  B0 = priors$B0
  c0 = priors$c0
  g0 = priors$g0
  G0 = priors$G0

  #empty objects to store parameters
  mu = matrix(NA_real_, nb_iter, K)
  sigma2 = matrix(NA_real_, nb_iter, K)
  eta = matrix(NA_real_, nb_iter, K)
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
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
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
    
    S = draw_allocations(pnorm)
    
    # 3. sample hyperparameters
    
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2[m, ]))
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))
    
    ## counter
    if(print){
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

#' @keywords internal
gibbs_SFM_poisson <- function(y,
                              K,
                              nb_iter,
                              priors = list(),
                              print = TRUE){
  
  # unpacking priors
  a0 = priors$a0
  A0 = priors$A0
  l0 = priors$l0
  L0 = priors$L0

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
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
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
      probs[,k] = eta[m,k]*dpois_unique(y,lambda[m,k])
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) 
    
    if (any(!is.finite(pnorm))) {
      stop("Non-finite count allocation probabilities; check support or numerical underflow.")
    }
    
    S = draw_allocations(pnorm)
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))
    
    ## counter
    if(print){
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

#' @keywords internal
gibbs_SFM_skew_n <- function(y,
                             K,
                             nb_iter,
                             priors = list(),
                             print = TRUE){
  
  # unpacking priors
  a0 = priors$a0
  A0 = priors$A0
  b0 = priors$b0
  c0 = priors$c0
  C0 = priors$C0
  g0 = priors$g0
  G0 = priors$G0
  D_xi = priors$D_xi
  D_psi = priors$D_psi
  
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
  e0 = a0/A0
  
  for (k in 1:K){
    xi[1, k] = mean(y[S[, k]==1])
  }
  sigma2 = rep(max(y)^2,K)
  beta = cbind(xi[1, ],0)
  
  zk = matrix(0,n_obs,1)
  cnt_update_e0 = 0
  
  for (m in 1:nb_iter){
    # set base priors
    ck = c0
    Ck = C0
    bk = b0
    Bk = solve(B0_inv)
    
    ## a.1 sample component proportion
    N = colSums(S)
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    probs = matrix(NA, n_obs, K)
    
    for (k in 1:K){
      
      if (N[k] > 0) {
        empty = FALSE
      } else {
        empty = TRUE
      }
      
      # sample z
      if(!empty){
        # allocate y
        yk = y[S[ ,k] == 1]
        
        # update z
        Ak = 1/(1 + beta[k, 2]^2/sigma2[k])
        ak = Ak*beta[k, 2]/sigma2[k]*(yk-beta[k, 1])
        zk <- rtnorm(N[k], ak, sqrt(Ak), lower=0) 
        
        Xk = matrix(c(rep(1, N[k]), zk), nrow=N[k])
      }
      
      ## a.1 sample Sigma
      if(!empty){
        eps = yk - Xk%*%beta[k, ]
        Ck = C0 + 0.5*sum(eps^2)
        ck = c0 + N[k]/2
      }
      sigma2[k] = 1/rgamma(1, ck, Ck)
      
      ## a.2 sample xi and psi jointly
      if(!empty){
        Bk = solve(crossprod(Xk)/sigma2[k] + B0_inv)
        bk = Bk%*%(crossprod(Xk, yk)/sigma2[k] + B0_inv%*%b0)
      }
      beta[k, ] = rmvnorm(1, bk, Bk)
      
      # storing
      xi[m, k] = beta[k, 1]
      omega[m, k] = sqrt(sigma2[k] + beta[k, 2]^2)
      alpha[m, k] = beta[k, 2]/sqrt(sigma2[k])
      probs[, k] = eta[m, k] * dsn(y, xi[m, k], omega[m, k], alpha[m, k])
    }
    
    # classification
    pnorm = probs/rowSums(probs)
    
    ## if the initial classification is bad then some data points won't be 
    # allocated to any components and some rows will be 
    # NAs (because if dividing by zero). We correct this by replacing NAs with
    # equal probabilities
    NA_id = which(is.na(pnorm[,1]))
    pnorm[NA_id, ] = 1/ncol(pnorm)
    
    S = draw_allocations(pnorm)
    
    # 3. sample hyperparameters
    ## a. sample C0
    C0 = rgamma(1, g0 + K*c0, G0 + sum(1/sigma2))
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))
    
    ## counter
    if(print){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc_result = cbind(eta, xi, omega, alpha, lp)
  colnames(mcmc_result) = 1:ncol(mcmc_result)
  
  for (i in 1:K){
    colnames(mcmc_result)[c(i, K+i, 2*K+i,3*K+i)] = c(paste0("eta", i),
                                                      paste0("xi", i),
                                                      paste0("omega", i),
                                                      paste0("alpha", i))
  }
  colnames(mcmc_result)[ncol(mcmc_result)] = "loglik"
  
  return(mcmc_result)
}

#' @keywords internal
gibbs_SFM_shifted_poisson <- function(y,
                         K,
                         nb_iter,
                         priors = list(),
                         print = TRUE){
  
  # unpacking priors
  a0 = priors$a0
  A0 = priors$A0
  l0 = priors$l0
  L0 = priors$L0
  
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
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    for (k in 1:K){
      
      yk = y[S[, k]==1]
      if (N[k] == 0) {
        kappa_m[k] = sample.int(max(y) + 1, 1) - 1
      } else {
        kapub = min(yk)
        pars_m = c(kappa_m[k], lambda_m[k])
        temp = draw_kap(yk, pars_m, kaplb = 0, kapub = kapub, dist = "shifted_poisson")
        kappa_m[k] = temp[[1]]
      }
      
      # Sample lambda from Gamma distribution 
      lambda_m[k] = rgamma(1, shape = sum(yk) - N[k]*kappa_m[k] + l0,
                           scale = 1/(N[k] + L0))
      
      # 
      probs[,k] = eta[m,k]*dpois_unique(y - kappa_m[k], lambda_m[k])
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs)
    
    if (any(!is.finite(pnorm))) {
      stop("Non-finite count allocation probabilities; check support or numerical underflow.")
    }
    
    S = draw_allocations(pnorm)
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))
    
    # storing
    lambda[m,] = lambda_m
    kappa[m, ] = kappa_m
    
    ## counter
    if(print){
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

# Gibbs sampler for the negative binomial
#' @keywords internal
gibbs_SFM_neg_binomial <- function(y,
                              K,
                              nb_iter,
                              nb_burnin,
                              priors = list(),
                              print = TRUE){
  
  # p in the code is the probability of failure
  # rtilde needs to be between 0 and 0.5
  # this implies r > 1

  # unpacking priors
  # priors for the negative binomial are implicit
  a0 = priors$a0
  A0 = priors$A0
  std_mh = rep(priors$sig, K)

  n_obs <- length(y)
  mu_max = max(y)
  if (mu_max <= 0) {
    stop("The negative binomial mean prior requires a positive full-data maximum.")
  }

  # storage matrices
  r_tilde = matrix(data=NA,nrow=nb_iter,ncol=K)
  mu = matrix(data=NA,nrow=nb_iter,ncol=K)  
  probs = matrix(data=NaN,nrow=n_obs,ncol=K) # Storage for probabilities  
  lp = matrix(0, nb_iter, 1) # Storage for log likelihood  
  eta = matrix(data = NA, nrow = nb_iter, ncol = K)
  
  # Initial conditions
  cl_y <- kmeans(y, centers = K, nstart = 30)
  batch_accept = rep(0, K)
  batch_proposals = rep(0, K)

  S <- matrix(0,length(y),K)
  
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  e0 = a0 / A0
  
  # initial value for the parameters
  mu_m <- cbind(t(cl_y$centers))
  mu_m[mu_m == 0] = mu_max / 2
  r_tilde_m <- rep(0.09, K)
  
  ## Sample lamda and S for each component, k=1,...,k
  for(m in 1:nb_iter){
    
    # Compute number of observations allocated in each component
    N = colSums(S)
    
    ## sample component proportion
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
    eta[m, ] = rdirichlet(1, e0 + N) 
    
    for (k in 1:K){
      
      yk = y[S[, k]==1]
     
      draw_pars = rneg_binom(yk, std_mh[k], mu_m[k], r_tilde_m[k], mu_max = mu_max)
      mu_m[k] = draw_pars[1]
      r_tilde_m[k] = draw_pars[2]
      if (m <= nb_burnin && N[k] > 0) {
        batch_accept[k] = batch_accept[k] + draw_pars[3]
        batch_proposals[k] = batch_proposals[k] + 1
      }
     # transform parameters
      r = (1-r_tilde_m[k])/r_tilde_m[k]
      p = mu_m[k] / (r + mu_m[k])
      
      #
      probs[,k]  = eta[m,k] * dnbinom_unique(y, r, 1 - p)
    }

    if (m <= nb_burnin && m %% 50 == 0) {
      std_mh = adapt_rtilde_proposals(std_mh, batch_accept, batch_proposals, m / 50)
      batch_accept[] = 0
      batch_proposals[] = 0
    }
    
    # 2. classification
    pnorm = probs/rowSums(probs) 
    
    if (any(!is.finite(pnorm))) {
      stop("Non-finite count allocation probabilities; check support or numerical underflow.")
    }
    
    S = draw_allocations(pnorm)
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))

    # storing
    mu[m,] = mu_m
    r_tilde[m, ] = r_tilde_m
    
    ## counter
    if(print){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }

  # output
  mcmc = cbind(eta, mu, r_tilde, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K) {
    colnames(mcmc)[c(i, K + i,  2*K+i)] = c(
      paste0("eta", i),
      paste0("mu", i),
      paste0("rtilde", i)
    )
  }
  
  colnames(mcmc)[ncol(mcmc)] = "loglik"

  # Return output   
  return(mcmc)
}

# Shifted negative binomial
#' @keywords internal
gibbs_SFM_shifted_neg_binomial <- function(y,
                         K,
                         nb_iter,
                         nb_burnin,
                         priors = list(),
                         print = TRUE){
  
  # unpacking priors
  a0 = priors$a0
  A0 = priors$A0
  std_mh = rep(priors$sig, K)
  
  n_obs <- length(y)
  mu_max = max(y)
  if (mu_max <= 0) {
    stop("The negative binomial mean prior requires a positive full-data maximum.")
  }
  
  # storage matrices
  kappa = matrix(data=NA,nrow=nb_iter,ncol=K) 
  r_tilde = matrix(data=NA,nrow=nb_iter,ncol=K)
  mu = matrix(data=NA,nrow=nb_iter,ncol=K)  
  probs = matrix(data=NaN,nrow=n_obs,ncol=K) # Storage for probabilities  
  lp = matrix(0, nb_iter, 1) # Storage for log likelihood  
  eta = matrix(data = NA, nrow = nb_iter, ncol = K)
  
  # Initial conditions
  cl_y <- kmeans(y, centers = K, nstart = 30)
  batch_accept = rep(0, K)
  batch_proposals = rep(0, K)

  S <- matrix(0,length(y),K)
  for (k in 1:K) {
    S[cl_y$cluster==k ,k] = 1
  }
  
  e0 = a0/A0
  
  # initial value for the parameters
  kappa_m = rep(0, K)
  mu_m <- cbind(t(cl_y$centers))
  mu_m[mu_m == 0] = mu_max / 2
  r_tilde_m <- rep(0.4, K)
  accept = rep(0, K)
  ## Sample lamda and S for each component, k=1,...,k
  for(m in 1:nb_iter){

    if (m <= nb_burnin){
      burnin = TRUE
    } else {
      burnin = FALSE
    }
    
    # Compute number of observations allocated in each component
    N = colSums(S)
    
    ## sample component proportion
    e0 = draw_e0(e0, a0, 1/A0, counts = N)[[1]]
    eta[m, ] = rdirichlet(1, e0 + N)

    for (k in 1:K){
      yk = y[S[, k] == 1]
      if (N[k]==0) {
        kappa_m[k] = sample.int(mu_max + 1, 1) - 1
      } else {
        draw_shift = draw_kap_mu(yk, kappa_m[k], mu_m[k], r_tilde_m[k], mu_max)
        kappa_m[k] = draw_shift$kappa
        mu_m[k] = draw_shift$mu
      }

      draw_pars = rneg_binom(yk - kappa_m[k], std_mh[k], mu_m[k], r_tilde_m[k], mu_max = mu_max)

      mu_m[k] = draw_pars[1]
      r_tilde_m[k] = draw_pars[2]
      mu_m[k] = draw_pars[1]
      r_tilde_m[k] = draw_pars[2]
      accept[k] = accept[k] + draw_pars[3]
      if (burnin && N[k] > 0) {
        batch_accept[k] = batch_accept[k] + draw_pars[3]
        batch_proposals[k] = batch_proposals[k] + 1
      }

      # transform parameters
      r = (1-r_tilde_m[k])/r_tilde_m[k]
      p = mu_m[k] / (r + mu_m[k])
      
      #
      probs[,k]  = eta[m,k] * dnbinom_unique(y - kappa_m[k], r, 1 - p)
      # Order the components based on kappa

    }

    if (burnin && m %% 50 == 0) {
      std_mh = adapt_rtilde_proposals(std_mh, batch_accept, batch_proposals, m / 50)
      batch_accept[] = 0
      batch_proposals[] = 0
    }

    # Ordering
    order_k = order(kappa_m)
    kappa_m = kappa_m[order_k]
    mu_m = mu_m[order_k]
    r_tilde_m = r_tilde_m[order_k]
    eta[m, ] = eta[m, order_k]
    probs = probs[, order_k, drop = FALSE]
    std_mh = std_mh[order_k]
    batch_accept = batch_accept[order_k]
    batch_proposals = batch_proposals[order_k]
    accept = accept[order_k]
    
    # 2. classification
    pnorm = probs/rowSums(probs)
    
    if (any(!is.finite(pnorm))) {
      stop("Non-finite count allocation probabilities; check support or numerical underflow.")
    }
    
    S = draw_allocations(pnorm)
    
    # compute log lik
    lp[m] = sum(log(rowSums(probs)))
    
    # storing
    mu[m,] = mu_m
    r_tilde[m, ] = r_tilde_m
    kappa[m, ] = kappa_m
    
    ## counter
    if(print){
      if(m %% (round(nb_iter / 10)) == 0){
        cat(paste(100 * m / nb_iter, ' % draws finished'), fill=TRUE)
      }
    }
  }
  
  # output
  mcmc = cbind(eta, kappa, mu, r_tilde, lp)
  colnames(mcmc) = 1:ncol(mcmc)
  
  for (i in 1:K) {
    colnames(mcmc)[c(i, K + i,  2*K+i, 3*K+i)] = c(
      paste0("eta", i),
      paste0("kappa", i),
      paste0("mu", i),
      paste0("rtilde", i)
    )
  }
  colnames(mcmc)[ncol(mcmc)] = "loglik"
  print(accept/nb_iter)
  # Return output   
  return(mcmc)
}

#' @keywords internal
check_priors <- function(priors, dist, data) {
  assert_that(all(is.finite(unlist(priors))),
              msg = "All priors should be finite.")
  
  # all
  priors$a0 = ifelse(is.null(priors$a0), 1, priors$a0)
  priors$A0 = ifelse(is.null(priors$A0), 200, priors$A0)
  
  assert_that(is.scalar(priors$a0), priors$a0 > 0, msg = "prior A0 should be a scalar")
  assert_that(is.scalar(priors$A0), priors$A0 > 0, msg = "prior A0 should be a positive scalar")
  
  
  if (dist == "shifted_poisson") {
    priors_labels = c("a0", "A0", "l0", "L0")
    
    priors$l0 = ifelse(is.null(priors$l0), 5, priors$l0)
    priors$L0 = ifelse(is.null(priors$L0), priors$l0 - 1, priors$L0)
  }
  
  if (dist == "poisson") {
    priors_labels = c("a0", "A0", "l0", "L0")

    priors$l0 = ifelse(is.null(priors$l0), 1.1, priors$l0)
    priors$L0 = ifelse(is.null(priors$L0), 1.1 / median(data), priors$L0)
  }
  
  if (dist %in% c("shifted_poisson", "poisson")) {
    assert_that(is.scalar(priors$L0), priors$L0 > 0, msg = "prior L0 should be a positive scalar")
    assert_that(is.scalar(priors$l0), priors$L0 > 0, msg = "prior l0 should be a positive scalar")
  }

  if (dist %in% c("shifted_neg_binomial", "neg_binomial")) {
    priors_labels = c("a0", "A0", "sig")

    priors$sig = ifelse(is.null(priors$sig), 0.1, priors$sig)
  }
  
  if (dist == "normal") {
    priors_labels = c("a0", "A0", "b0", "B0", "c0", "g0", "G0")
    
    priors$b0 = ifelse(is.null(priors$b0), median(data), priors$b0)
    priors$B0 = ifelse(is.null(priors$B0), (max(data) - min(data))^2, priors$B0)
    priors$c0 = ifelse(is.null(priors$c0), 2.5, priors$c0)
    priors$g0 = ifelse(is.null(priors$g0), 0.5, priors$g0)
    priors$G0 = ifelse(is.null(priors$G0), 100*priors$g0/priors$c0/priors$B0, priors$G0)
    
    assert_that(is.scalar(priors$b0), msg = "prior b0 should be a scalar")
    assert_that(is.scalar(priors$B0), priors$B0 > 0, msg = "prior B0 should be a positive scalar")
    assert_that(is.scalar(priors$c0), priors$c0 > 0, msg = "prior c0 should be a positive scalar")
    assert_that(is.scalar(priors$g0), priors$g0 > 0, msg = "prior g0 should be a positive scalar")
    assert_that(is.scalar(priors$G0), priors$G0 > 0, msg = "prior G0 should be a positive scalar")
  }
  
  if (dist == "skew_normal") {
    priors_labels = c("a0", "A0", "b0", "c0", "C0", "g0", "G0", "D_xi", "D_psi")
    
    priors$b0 = ifelse(is.null(priors$b0), median(data), priors$b0)
    priors$c0 = ifelse(is.null(priors$c0), 2.5, priors$c0)
    priors$C0 = ifelse(is.null(priors$C0), 0.5*var(data), priors$C0)
    priors$g0 = ifelse(is.null(priors$g0), 0.5, priors$g0)
    priors$G0 = ifelse(is.null(priors$G0), priors$g0/(0.5*var(data)), priors$G0)
    priors$D_xi = ifelse(is.null(priors$D_xi), 1, priors$D_xi)
    priors$D_psi = ifelse(is.null(priors$D_psi), 1, priors$D_psi)
    
    assert_that(is.scalar(priors$b0), msg = "prior b0 should be a scalar")
    assert_that(is.scalar(priors$D_xi), priors$D_xi > 0, msg = "prior D_xi should be a positive scalar")
    assert_that(is.scalar(priors$D_psi), priors$D_psi > 0, msg = "prior D_psi should be a positive scalar")
    assert_that(is.scalar(priors$c0), priors$c0 > 0, msg = "prior c0 should be a positive scalar")
    assert_that(is.scalar(priors$C0), priors$C0 > 0, msg = "prior C0 should be a positive scalar")
    assert_that(is.scalar(priors$g0), priors$g0 > 0, msg = "prior g0 should be a positive scalar")
    assert_that(is.scalar(priors$G0), priors$G0 > 0, msg = "prior G0 should be a positive scalar")
  }
  
  isnot = which(!names(priors) %in% priors_labels)
  
  if (length(isnot) > 0) {
    warning(paste("prior(s)", names(priors)[isnot], "not needed when estimating a", dist, "mixture."))
  }
  
  return (priors[priors_labels])
}

## functions used in the SFM MCMC algorithm

draw_allocations <- function(probabilities) {
  if (!is.matrix(probabilities) || ncol(probabilities) == 0L ||
      any(!is.finite(probabilities)) || any(probabilities < 0)) {
    stop("Allocation probabilities must be a finite, non-negative matrix with at least one component.")
  }
  totals = rowSums(probabilities)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("Allocation probabilities must have positive finite row sums.")
  }
  count = nrow(probabilities)
  components = ncol(probabilities)
  allocations = matrix(0L, nrow = count, ncol = components)
  if (count == 0L) return(allocations)
  if (components == 1L) {
    allocations[, 1] = 1L
    return(allocations)
  }
  thresholds = runif(count) * totals
  cumulative = numeric(count)
  selected = integer(count)
  last_positive = integer(count)
  for (component in seq_len(components)) {
    last_positive[probabilities[, component] > 0] = component
    cumulative = cumulative + probabilities[, component]
    selected[selected == 0L & thresholds < cumulative] = component
  }
  unresolved = selected == 0L
  selected[unresolved] = last_positive[unresolved]
  allocations[cbind(seq_len(count), selected)] = 1L
  allocations
}

dpois_unique <- function(x, lambda, log = FALSE) {
  if (length(lambda) != 1L || length(log) != 1L) return(dpois(x, lambda, log = log))
  values = unique(x)
  probabilities = dpois(values, lambda, log = log)
  probabilities[match(x, values)]
}

dnbinom_unique <- function(x, ...) {
  if (any(lengths(list(...)) != 1L)) return(dnbinom(x, ...))
  values = unique(x)
  probabilities = dnbinom(values, ...)
  probabilities[match(x, values)]
}

adapt_rtilde_proposals <- function(std, accepted, proposed, batch) {
  acceptance_rate = accepted / pmax(proposed, 1)
  direction = as.numeric(proposed > 0 & acceptance_rate > 0.4) -
    as.numeric(proposed > 0 & acceptance_rate < 0.3)
  std * exp(min(0.1, 1 / sqrt(batch)) * direction)
}

#' @keywords internal
rneg_binom <- function(Y, std, d_mu, d_rtilde, mu_max) {

  #Y = Y[Y >= 0]
  
  # Inputs
  # Y: Nx1 vector of observations
  # std: scalar valued std dev for MH step 
  # Initial values of:
  #   d_mu: scalar valued mean
  #   d_rtilde: scalar valued rtilde
  
  # Output
  # Updated values of:
  #   d_mu: scalar valued mean
  #   d_rtilde: scalar valued rtilde
  
  # Sample
  # Define useful terms
  N <- length(Y) # Sample size
  if (length(mu_max) != 1 || !is.finite(mu_max) || mu_max <= 0) {
    stop("mu_max must be a positive finite scalar.")
  }
  if (N == 0) {
    return(c(runif(1, 0, mu_max), runif(1, 0, 0.5), 0, std))
  }
  
  # Define NB terms
  d_r <- (1 - d_rtilde) / d_rtilde
  d_p <- d_mu / (d_r + d_mu)
  
  # Simulate vector of lambda_i (i=1,2,....,n) from
  # Gamma(alpha = y_i + r, beta = 1 / p) distribution
  d_b <- d_p
  v_lambda <- rgamma(N, shape = Y + d_r, scale = d_b)
  
  # Simulate 1 / mu from truncated version on interval (1 / M, infinity) of
  # Gamma(alpha = rn - 1, beta = r * sum(lambda_i))
  d_a <- d_r * N - 1
  d_b <- 1 / (d_r * sum(v_lambda))
  d_log_tail_at_1_over_m <- pgamma(1 / mu_max, shape = d_a, scale = d_b,
                                  lower.tail = FALSE, log.p = TRUE)

  d_log_u <- d_log_tail_at_1_over_m + log(runif(1))
  d_mu <- 1 / qgamma(d_log_u, shape = d_a, scale = d_b,
                     lower.tail = FALSE, log.p = TRUE)
  
  if (d_mu == 0) d_mu <- 1e-8
  
  # Update value of d_p and compute loglikelihood (excluding the term 
  # -sum log(y_i!) that only depends on the data) of current draw (mu, r~)
  d_p <- d_mu / (d_r + d_mu)
  #d_log_l <- d_r * N * log(1 - d_p) - N * lgamma(d_r) + sum(lgamma(Y + d_r) + Y * log(d_p))
  
  d_log_l = sum(dnbinom_unique(Y, size = d_r, prob = 1 - d_p, log = TRUE))

  d_rtilde_candidate = rnorm(1, mean = d_rtilde, sd = std)
  accept = 0
  if (d_rtilde_candidate > 0 & d_rtilde_candidate < 0.5) {
    d_r_candidate <- (1 - d_rtilde_candidate) / d_rtilde_candidate
    d_p_candidate <- d_mu / (d_r_candidate + d_mu)
    #d_log_l_candidate <- d_r_candidate * N * log(1 - d_p_candidate) -
    #  N * lgamma(d_r_candidate) + sum(lgamma(Y + d_r_candidate) + Y * log(d_p_candidate))
    d_log_l_candidate = sum(dnbinom_unique(Y, size = d_r_candidate, prob = 1 - d_p_candidate, log = TRUE))

    d_log_acceptance_probability <- min(d_log_l_candidate - d_log_l, 0)
    
    d_u_rw_mh <- runif(1)
    if (log(d_u_rw_mh) < d_log_acceptance_probability) {
      accept = 1
      d_r <- d_r_candidate
      d_rtilde <- d_rtilde_candidate
      d_p <- d_p_candidate
    }
  }
  
  return(c(d_mu, d_rtilde, accept, std))
}

draw_kap_mu <- function(x, kappa, mu, r_tilde, mu_max) {
  component_mean = kappa + mu
  shifts = seq.int(0, min(x, mu_max))
  means = component_mean - shifts
  supported = means > 0 & means <= mu_max
  shifts = shifts[supported]
  means = means[supported]
  if (!length(shifts)) {
    stop("Joint shift-mean update has no candidates within the conditional support.")
  }

  shape = (1 - r_tilde) / r_tilde
  log_probabilities = vapply(seq_along(shifts), function(candidate) {
    sum(dnbinom_unique(x - shifts[candidate], size = shape, mu = means[candidate], log = TRUE))
  }, numeric(1))
  maximum = max(log_probabilities)
  if (!is.finite(maximum)) {
    stop("Non-finite joint shift-mean conditional probabilities.")
  }
  selected = if (length(shifts) == 1) 1L else {
    sample.int(length(shifts), 1, prob = exp(log_probabilities - maximum))
  }
  list(kappa = shifts[selected], mu = means[selected])
}

# Posterior of kappa 
#' @keywords internal
post_kap <- function(x, pars, dist){

  kappa = pars[1]

  if (dist == "shifted_poisson") {
    lambda = pars[2]
    result = sum(dpois_unique(x - kappa, lambda = lambda, log = TRUE))
  }

  if (dist == "shifted_neg_binomial") {
    r_tilde = pars[2]
    mu = pars[3]
    
    # transform parameters
    r = (1-r_tilde)/r_tilde
    p = mu / (r + mu)
    
    result = sum(dnbinom(x - kappa, r, 1 - p, log=TRUE))
    #result = neg_binom_ll(x - kappa, r, p)
    #result = exp(sum(lgamma(x - kappa + r) - lfactorial(x - kappa)) + sum(x - kappa) * log(p))
  }

  return (result)
}

# Draw kappa from posterior using MH step
#' @keywords internal
draw_kap <- function(x,pars,kaplb,kapub, dist) {
  KAPPA = pars[1]
  if (KAPPA < kaplb || KAPPA > kapub) {
    stop("Current shift is outside its conditional support.")
  }
  KAPPAhat = KAPPA + sample(c(-1,0,1), 1)
  if (KAPPAhat < kaplb || KAPPAhat > kapub) {
    return(list(KAPPA, 0))
  }
  pars_hat = c(KAPPAhat, pars[-1])

  log_ratio = post_kap(x, pars_hat, dist) - post_kap(x, pars, dist)
  if (!is.finite(log_ratio)) stop("Non-finite shift log acceptance ratio.")
  rand = runif(1, min = 0, max = 1) # Random draw from unifrom (0,1)
  if (log(rand) < min(0, log_ratio)) {
    KAPPA = KAPPAhat # update
    acc = 1 # indicate update
  } else {
    # don't update
    acc = 0 # indicate failure to update
  }
  
  out <- list(KAPPA, acc) # Store output in list
  return(out)           # Return output
}

# Posterior of e0 - Unnormalized target pdf
#' @keywords internal
post_e0 <- function(e0,nu0_p,S0_p,p = NULL, log = FALSE, counts = NULL) {
  if (!is.finite(e0) || e0 <= 0) {
    return(if (log) -Inf else 0)
  }
  if (!is.null(counts)) {
    if (!length(counts) || any(!is.finite(counts)) || any(counts < 0 | counts %% 1 != 0)) {
      stop("Allocation counts must be non-negative finite integers.")
    }
    result = dgamma(e0, shape = nu0_p, scale = S0_p, log = TRUE) +
      lgamma(length(counts) * e0) - lgamma(sum(counts) + length(counts) * e0) +
      sum(lgamma(counts + e0) - lgamma(e0))
    return(if (log) result else exp(result))
  }
  if (any(!is.finite(p)) || any(p <= 0)) {
    stop("Concentration updates require strictly positive finite mixture weights.")
  }
  K = length(p) # Number of components
  result = dgamma(e0, shape = nu0_p, scale = S0_p, log = TRUE) +
    lgamma(K * e0) - K * lgamma(e0) + (e0 - 1) * sum(log(p))
  return(if (log) result else exp(result))
}

# Draw from the unnormalized target pdf for hyperparameter e0
#' @keywords internal
draw_e0 <- function(e0,nu0,S0,p = NULL, counts = NULL){
  # Define terms
  e0hat = e0 + rnorm(1, 0, 0.1) # Draw a candidate from Random walk proposal
  if (!is.finite(e0hat) || e0hat <= 0) {
    return(list(e0, 0))
  }
  log_ratio = post_e0(e0hat, nu0, S0, p, log = TRUE, counts = counts) -
    post_e0(e0, nu0, S0, p, log = TRUE, counts = counts)
  if (!is.finite(log_ratio)) stop("Non-finite concentration log acceptance ratio.")
  # MH Step
  rand = runif(1, min = 0, max = 1) # Random draw from unifrom (0,1)
  if (log(rand) < min(0, log_ratio)) {
    e0 = e0hat; # update
    acc = 1; # indicate update
  } else{
    e0 = e0; # don't update
    acc = 0;  # indicate failure to update
  }
  out <- list(e0, acc) # Store output in list
  return(out)           # Return output
}

# for the negative binomial
neg_binom_ll <- function(y, r, p) {
  # Evaluates the negative binomial log-likelihood using the paramterization
  # where p is the probability of success, and r is the number of successes,
  # i.e.,
  # P(Y=y)=(gamma(y + r)/(gamma(y+1)*gamma(r)))*p^r*(1-p)^y
  #
  # Inputs:
  #   y     - Nx1 vector of observed counts (non-negative integers).
  #   p     - probability of success (real number in (0,1)).
  #   r     - Dispersion parameter (positive real number).
  #
  # Output:
  #   logLikelihood - The negative binomial log-likelihood.

  # Check if any adjusted counts are negative
  if (any(y < 0)) {
    logLikelihood <- -Inf # Return -Inf to invalidate such observations
  } else {
    # Compute the log-likelihood
    term1 <- lgamma(y + r) - lgamma(r) - lgamma(y + 1)
    term2 <- r * log(p)
    term3 <- y * log(1 - p)

    logLikelihood <- sum(term1 + term2 + term3)

  }

  return(logLikelihood)
}