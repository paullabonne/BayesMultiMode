% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_mode.R
\name{bayes_mode}
\alias{bayes_mode}
\title{Mode inference using post-processed SFM MCMC draws.}
\usage{
bayes_mode(theta_draws, y, mixt = "shifted_poisson")
}
\arguments{
\item{theta_draws}{a (M x 3xJb) matrix. Output of `sfm_mcmc_spmix()` giving MCMC parameter draws after burn-in and discarding empty components.}

\item{y}{(a vector of integers) Observations used to fit the model.}

\item{mixt}{(a string) giving the mixture distribution. Default is "shifted_poisson".}
}
\value{
A list containing:
\itemize{
  \item Prob_unimod : Posterior probability of unimodality. (1-Prob_unimod) is equal to the posterior probability of multimodality.
  \item table_nb_modes : Possible number of modes and posterior probability for each of those.
  \item table_locations : Possible locations of modes and posterior probability for each of those.
  \item A list of graphs showing : 
  \itemize{
       \item 1: The posterior probability of multimodality;
       \item 2: Possible number of modes and posterior probability for each of those;
       \item 3: Possible locations of modes and posterior probability for each of those;
}
}
}
\description{
Computes the number of modes, their locations and posterior probabilities.
}
\examples{
# Example with simulated data ================================================
#set seed for random number generation
set.seed(1) 

# Set the parameters for drawing from a two-component shifted Poisson mixture
p1 = 0.3
p2 = 1-p1
kap1 = 3
kap2 = 0
lam1 = 1
lam2 = 0.5
length_data = 70

# Generate data
y <- c(rpois(length_data*p1, lam1)+kap1, rpois(length_data*p2, lam2)+kap2)

# Set parameters for the SFM MCMC estimation
M = 1000 # Number of MCMC iterations 
Jmax = 4 # Maximum number of mixture components

# Estimation with SFM MCMC
sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)

# Proportion of draws burned in
S = 0.5

# Post processing
post_sfm_mcmc = post_sfm_mcmc(sfm_mcmc, S=S)

# Mode inference
bayes_mode(post_sfm_mcmc$theta_draws_slim,y)

# Example with DNA data =====================================================
\donttest{
y = d4z4
M = 5000 # Number of MCMC iterations 
Jmax = 10 # Maximum number of mixture components

# Estimation with SFM MCMC

sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
# Proportion of draws burned in
S = 0.5

# Mode inference
bayes_mode(post_sfm_mcmc$theta_draws_slim,y)
}
}
\references{
\insertRef{basturk_bayes_2021}{BayesMultiMode}
}
