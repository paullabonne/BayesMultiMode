#' Estimation of a mixture of shifted Poisson distributions.
#' 
#' Bayesian estimation of a mixture of shifted Poisson distributions using a Sparse Finite Mixture MCMC algorithm.
#' @param y (a vector of integers) Observations used to fit the model.
#' @param Jmax (an integer) Maximum number of mixture components.
#' @param M (an integer) Number of MCMC iterations.
#' @param prt print intermediate of the MCMC estimation ? default = TRUE.
#' @returns 
#' A list containing:
#' \itemize{
#'   \item theta_draws : a (M x 3xJmax) matrix. Parameter draws from the posterior distribution at each MCMC iteration.
#'   \item p_draws : a (M x Jmax) matrix. Posterior draw of mixture weights at each MCMC iteration.
#'   \item kappa_draws : a (M x Jmax) matrix. Posterior draws of kappa at each MCMC iteration.
#'   \item lambda_draws : a (M x Jmax) matrix. Posterior draws of lambda at each MCMC iteration.
#'   \item snj : a (M x Jmax) matrix. Number of observations allocated to each components at each MCMC iteration.
#'   \item se0 : a vector of size M. Concentration parameter from (symmetric) Dirichlet distribution at each MCMC iteration.
#'   \item dist : a string indicating the distribution used in the mixture.
#'   \item y, Jmax and M, given as input.
#' }
#' 
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{malsiner-walli_model-based_2016}{BayesMultiMode}\cr
#' \insertRef{basturk_bayes_2021}{BayesMultiMode}
#' @examples
#'# Example with simulated data ================================================
#' #set seed for random number generation
#' set.seed(1) 
#' 
#' # Set the parameters for drawing from a two-component shifted Poisson mixture
#' p1 = 0.3
#' p2 = 1-p1
#' kap1 = 3
#' kap2 = 0
#' lam1 = 1
#' lam2 = 0.5
#' length_data = 70
#' 
#' # Generate data
#' y <- c(rpois(length_data*p1, lam1)+kap1, rpois(length_data*p2, lam2)+kap2)
#' 
#' # Set parameters for the SFM MCMC estimation
#' M = 1000 # Number of MCMC iterations
#' Jmax = 4 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' 
#' # Example with DNA data =====================================================
#' \donttest{
#' y = d4z4
#' M = 5000 # Number of MCMC iterations 
#' Jmax = 10 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' }
#' @importFrom gtools rdirichlet
#' @importFrom stats density dgamma dpois rgamma rmultinom rnorm runif
#' 
#' @export

sfm_mcmc_spmix <- function(y, Jmax, M, prt = TRUE){
  # initialization 
  tmp <- runif(Jmax)
  theta.init <- c(tmp[-Jmax] / sum(tmp), rep(min(y),Jmax), rep(1,Jmax))
  names(theta.init) = c(paste("prob ",1:(Jmax-1),sep=""),
                        paste("kappa ",1:Jmax,sep=""),
                        paste("lambda ",1:Jmax,sep=""))
  theta0 <- theta.init
  
  # Error Messages  
  if(round(Jmax) != Jmax | Jmax < 1){
    stop("max number of mixture components should be integer >= 1")
  }
  if(!is.vector(y)){
    stop("data 'y' should be a vector")
  }

  # initial parameters
  ind = Jmax
  p0 <- theta0[1:(Jmax-1)];
  kappa0   <- theta0[ind:(ind+Jmax-1)]; ind=ind+Jmax
  lambda0 <- theta0[ind:(ind+Jmax-1)]; 
  p <- c(p0,(1-sum(p0))); lam <- lambda0; kap <- kappa0
  N <- length(y)     
  lambda.rng <- kappa.rng <-c(0,max(y))
  
  if(any(lambda0<=0))
    stop("'lambda0' should be > 0 for all clusters");
  if(any(kappa0>min(y)))
    stop("'kappa0' should be <= min(y) for all clusters");
  
  ## Priors
  # kappa_j - uniform prior where upper bound varies in MCMC
  kaplb = 0 # Lower bound for kappa
  
  # lambda_j - independent gamma prior - symmetric prior across all components
  nu0 = 5
  invS0 = (nu0-1)
  
  # e0 - independent gamma prior: Eq (3) in MW, FS and Grun
  nu0_p = 1 
  S0_p = 1/200
  
  # Initial conditions
  lam = matrix(data=lam,nrow=1,ncol=Jmax) # lambda 
  kap = matrix(data=kap,nrow=1,ncol=Jmax) # kappa 
  p = matrix(data=p,nrow=1,ncol=Jmax) # probabilities
  Si = t(rmultinom(n = N,size=1,prob=p)) # component indicator
  alpha0 = matrix(data=nu0_p*S0_p,nrow=1,ncol=Jmax) # concentration parameter - note mean of gamma distribution = nu0_p*S0_p
  
  # storage matrices
  slam = matrix(data=NA,nrow=M,ncol=Jmax) # lambda
  skap = matrix(data=NA,nrow=M,ncol=Jmax) # kappa
  sP = matrix(data=NA,nrow=M,ncol=Jmax)   # probabilities
  snj = matrix(data=NA,nrow=M,ncol=Jmax)  # number of draws in each component
  se0 = matrix(data=NA,nrow=M,ncol=1) # Symmetric Dirichlet hyperparameter
  sJ = matrix(data=NA,nrow=M,ncol=1) # Number of non-empty components
  
  # Counters for MH Steps
  cnt_update_kap = matrix(data=0,nrow=1,ncol=Jmax)  # counter for MH step update of kappa
  cnt_update_e0 = 0 # counter for MH step update of concentration parameter
  
  for(m in 1:M){
    ## Sample kappa, lamda and S for each component, j=1,...,k
    probs = matrix(data=NaN,nrow=N,ncol=Jmax) # Storage for probabilities  
    for (j in 1:Jmax){
      # Set-up
      indj = (Si[,j]==1) # Find all S_i = j
      nj = sum(indj)
      if (nj==0) {
        yj = 0
      } else {
        yj = y[indj]
      }
      
      
      # Sample kappa using MH Step
      kapub = min(yj)
      if (length(kapub) == 0) {# Set to zero if component is empty
        kap[,j] = 0; 
      } else if (kapub < kap[,j]) {# Set to upper bound if outside boudary
        kap[,j] = kapub
      } else {
        temp = draw_kap(yj,lam[,j],kap[,j],kaplb,kapub) # Draw kappa from MH step
        kap[,j] = as.numeric(temp[1])
        cnt_update_kap[,j] = cnt_update_kap[,j] + as.numeric(temp[2])
      }
      
      # Sample lambda from Gamma distribution - change if prior is non-Gamma
      nup = (sum(yj)-nj*kap[,j]) + nu0
      Sp = 1/(nj + invS0)
      lam[,j] = rgamma(1, shape = nup, scale = Sp)
      
      # Sample component indicator from multinomial
      probs[,j] = p[,j]*spoisspdf(y,lam[,j],kap[,j])        
    }
    
    pnorm = probs/replicate(Jmax,rowSums(probs))
    for (jj in 1:N){
      Si[jj,] = t(rmultinom(n = 1,size=1,prob=pnorm[jj,])) # component indicator
    }
    
    ## Sample component probabilities from Dirichlet
    nj2 = colSums(Si) # number of obs in each component
    alphastar = alpha0 + nj2 # posterior concentration parameter
    p = rdirichlet(1,alphastar) 
    
    ## Sample component probabilities hyperparameters: alpha0, using RWMH step  
    e0 = alpha0[1] # Symmetric parameter so choose 1 as default
    temp = draw_e0(e0,nu0_p,S0_p,p)
    e0 = as.numeric(temp[1])
    ind_update_e0 = as.numeric(temp[2])
    alpha0 = e0*matrix(data=1,nrow=1,ncol=Jmax) # probabilities
    cnt_update_e0 = cnt_update_e0 + ind_update_e0
    
    ## Ordering based on the expected value (mu = kappa+lambda)
    mu = kap+lam
    matrix_temp = t(rbind(mu,p,kap,lam,nj2))
    matrix_temp = matrix_temp[order(matrix_temp[,1]),]

    ## Store results  
    slam[m,] = matrix_temp[,4]    
    skap[m,] = matrix_temp[,3]     
    sP[m,] = matrix_temp[,2] 
    snj[m,] = matrix_temp[,5] 
    se0[m,] = e0  
    sJ[m,] = Jmax - sum(is.na(slam[m,]))
    
    # Set number of data points in empty components to NaN value
    sP[m, snj[m,]==0] = NA 
    skap[m, snj[m,]==0] = NA 
    slam[m, snj[m,]==0] = NA
    
    colnames(sP) = 1:ncol(sP)
    colnames(skap) = 1:ncol(skap)
    colnames(slam) = 1:ncol(slam)
    
    ## counter
    if(prt){
      if(m %% (round(M / 10)) == 0){
        naccept_e0 = 100*cnt_update_e0/m
        cat(paste(100 * m / M, ' % draws finished. Accept. prob of e0 =', round(naccept_e0), 'percent'), fill=TRUE)
      }
    }
  }
  
  # Store output
  theta_draws = cbind(sP, skap, slam)
  
  # Return output   
  return(list(theta_draws = theta_draws,
              p_draws = sP, kappa_draws = skap, lambda_draws = slam,
              snj=snj, se0 = se0,
              mixt = "shifted_poisson",
              y = y ,Jmax = Jmax, M = M))
}
