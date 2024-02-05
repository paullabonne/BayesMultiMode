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
    mcmc = gibbs_SFM_sp(y, K, nb_iter, priors, print)
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
gibbs_SFM_sp <- function(y,
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
    priors$L0 = ifelse(is.null(priors$L0), 1.1/median(data), priors$L0)
  }
  
  if (dist %in% c("shifted_poisson", "poisson")) {
    assert_that(is.scalar(priors$L0), priors$L0 > 0, msg = "prior L0 should be a positive scalar")
    assert_that(is.scalar(priors$l0), priors$L0 > 0, msg = "prior l0 should be a positive scalar")
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

# Posterior of kappa
#' @keywords internal
post_kap <- function(x,LAMBDA,KAPPA) {
  n = length(x) # Number of elements in the component
  result <-  exp(-LAMBDA)*(LAMBDA^(sum(x)-n*KAPPA))/prod(factorial(x-KAPPA)) 
}

# Draw kappa from posterior using MH step
#' @keywords internal
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

# Posterior of e0 - Unnormalized target pdf
#' @keywords internal
post_e0 <- function(e0,nu0_p,S0_p,p) {
  K = length(p) # Number of components
  result <-  dgamma(e0,shape = nu0_p, scale = S0_p)*gamma(K*e0)/(gamma(e0)^K)*((prod(p))^(e0-1))
}

# Draw from the unnormalized target pdf for hyperparameter e0
#' @keywords internal
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