# BayesMultiMode
An R package for testing multimodality in discrete distributions using Bayesian methods

### Installing the development version of BayesMultiMode
```{r}
#install.packages("devtools") #if needed
library(devtools)

devtools::install_github("paullabonne/BayesMultiMode")
```

### Generating data
```{r}
set.seed(1)
p1 = 0.3
p2 = 1-p1
kap1 = 3
kap2 = 0
lam1 = 1
lam2 = 0.5
length_data = 70
simulated_data <- c(rpois(length_data*p1, lam1)+kap1, rpois(length_data*p2, lam2)+kap2)
```

### Running the estimation on simulated or DNA data
```{r}
# Select DNA data :
y = d4z4

# Or select simulated data :
# y = simulated_data
```

### Setting parameters for MCMC
```{r}
# Number of MCMC iterations 
M = 5000 

# Proportion of draws to discard as burnin
S = 0.5 

# Maximum number of mixture components 
Jmax = 10
```

### Estimation with SFM MCMC
```{r}
#Bayesian estimation
sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)

#Plots of the estimation output
plots_mcmc(sfm_mcmc,S)

# Simple post-processing
post_sfmmcmc = post_sfm_mcmc(sfm_mcmc,S)
```

### Inference on modes
```{r}
sfm_mcmc_modes = bayes_mode(post_sfmmcmc$theta_draws_slim,y)
sfm_mcmc_modes$graphs
```
