Compatibility with other R packages
================

``` r
library(BayesMultiMode)
library(posterior)
```

## Bayesian estimation and mode inference

In the examples presented below external R packages are used for
Bayesian estimation of mixture models while `BayesMultiMode` is used for
mode inference.

### rjags

#### Estimation

``` r
set.seed(123)

library(rjags)

# Assuming you have your data in 'data_vector'
# Specify the model
model_string = "
  model {
    for (i in 1:N) {
      data_vector[i] ~ dnorm(mu[component[i]], tau[component[i]])
      component[i] ~ dcat(theta[])
    }

    for (j in 1:K) {
      mu[j] ~ dnorm(0, 0.01)
      tau[j] ~ dgamma(1, 1)
    }

    theta ~ ddirch(rep(1, K))
  }
"

K = 2
# Data for JAGS
y = c(rnorm(100,0,1),
      rnorm(100,5,1))

data_jags = list(data_vector = y, N = length(y), K = K)

# Initialize parameters
inits = function() {
  list(mu = rnorm(K, 0, 10),
       tau = rgamma(K, 1, 1),
       theta = rep(1 / K, K))
}

# Run the model
model = jags.model(textConnection(model_string), data = data_jags, inits = inits)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 200
    ##    Unobserved stochastic nodes: 205
    ##    Total graph size: 811
    ## 
    ## Initializing model

``` r
update(model, 1000)  # Burn-in
fit = coda.samples(model, variable.names = c("mu", "tau", "theta"), n.iter = 2000)
```

#### Create a BayesMixture object

``` r
fit_mat = as_draws_matrix(fit)
bmix = bayes_mixture(mcmc = fit_mat,
                        data = y,
                        burnin = 0, # the burnin has already been discarded
                        dist = "normal",
                        vars_to_keep = c("theta", "mu", "tau"),
                        vars_to_rename = c("eta" = "theta",
                                           "sigma" = "tau"))

# plot the mixture
plot(bmix)
```

![](external_comp_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# mode estimation
bayesmode = bayes_mode(bmix)

# plot mode inference
plot(bayesmode)
```

![](external_comp_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# Summary of mode inference
summary(bayesmode)
```

    ## Posterior probability of multimodality is 1 
    ## 
    ## Snapshot of inference results on the number of modes:
    ##   tb_nb_modes (matrix, dim 1x2): 
    ##      number of modes posterior probability
    ## [1,]               2                     1
    ## 
    ## Snapshot of inference results on mode locations:
    ##   table_location (matrix, dim 55x2): 
    ##      mode location posterior probability
    ## [1,]          -0.2                0.0125
    ## [2,]          -0.1                0.0875
    ## [3,]           0.0                0.3000
    ## [4,]           0.1                0.0000
    ## [5,]           0.2                0.1915
    ## [6,]           0.3                0.0325
    ## ... (49 more rows)

### rstan

#### Estimation

``` r
set.seed(123)
library(rstan)

normal_mixture_model <- "
data {
  int<lower=0> N;         // number of data points
  real y[N];              // observed data
  int<lower=1> K;         // number of mixture components
}

parameters {
  simplex[K] theta;       // mixing proportions
  ordered[K] mu;          // means of the Gaussian components
  vector<lower=0>[K] sigma; // standard deviations of the components
}

model {
  vector[K] log_theta = log(theta); // cache log calculation
  
  // Priors
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);

  // Likelihood
  for (n in 1:N) {
    vector[K] log_lik;
    for (k in 1:K) {
      log_lik[k] = log_theta[k] + normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(log_lik);
  }
}
"

y = c(rnorm(100,0,1),
      rnorm(100,5,1))

data_list <- list(
  N = length(y),
  y = y,
  K = 2  # Assuming a two-component mixture
)

fit <- stan(model_code = normal_mixture_model, data = data_list, iter = 2000, chains = 1)
```

    ## 
    ## SAMPLING FOR MODEL '07e8533f6a3188d2d50eb989867b63d7' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 4.9e-05 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.49 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 0.25367 seconds (Warm-up)
    ## Chain 1:                0.194892 seconds (Sampling)
    ## Chain 1:                0.448562 seconds (Total)
    ## Chain 1:

#### Create a BayesMixture object

``` r
fit_mat = as_draws_matrix(fit)
bmix = bayes_mixture(mcmc = fit_mat,
                        data = y,
                        burnin = 0, # the burnin has already been discarded
                        dist = "normal",
                        vars_to_keep = c("theta", "mu", "sigma"),
                        vars_to_rename = c("eta" = "theta"))

# plot the mixture
plot(bmix)
```

![](external_comp_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# mode estimation
bayesmode = bayes_mode(bmix)

# plot mode inference
plot(bayesmode)
```

![](external_comp_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# Summary of mode inference
summary(bayesmode)
```

    ## Posterior probability of multimodality is 1 
    ## 
    ## Snapshot of inference results on the number of modes:
    ##   tb_nb_modes (matrix, dim 1x2): 
    ##      number of modes posterior probability
    ## [1,]               2                     1
    ## 
    ## Snapshot of inference results on mode locations:
    ##   table_location (matrix, dim 55x2): 
    ##      mode location posterior probability
    ## [1,]          -0.2                 0.007
    ## [2,]          -0.1                 0.097
    ## [3,]           0.0                 0.308
    ## [4,]           0.1                 0.000
    ## [5,]           0.2                 0.194
    ## [6,]           0.3                 0.031
    ## ... (49 more rows)

### bayesmix

#### Estimation

``` r
set.seed(123)
library(bayesmix)

## taken from the examples of JAGSrun()
data("fish", package = "bayesmix")
y = unlist(fish)
prefix <- "fish"
variables <- c("mu","tau","eta")
k <- 3
modelFish <- BMMmodel(k = k, priors = list(kind = "independence",
                                           parameter = "priorsFish", hierarchical = "tau"))
controlFish <- JAGScontrol(variables = c(variables, "S"),
                           n.iter = 2000, burn.in = 1000)
z1 <- JAGSrun(fish, prefix, model = modelFish, initialValues = list(S0 = 2),
              control = controlFish, cleanup = TRUE, tmp = FALSE)
```

    ## Compiling model graph
    ##    Declaring variables
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 256
    ##    Unobserved stochastic nodes: 264
    ##    Total graph size: 1044
    ## 
    ## Initializing model

``` r
fit = z1$results
```

#### Create a BayesMixture object

``` r
# JAGSrun seems to report the variance (sigma2) rather the standard deviation so we need a new pdf :
pdf_func <- function(x, pars) {
  dnorm(x, pars["mu"], sqrt(pars["sigma"]))
}

bmix = bayes_mixture(mcmc = fit,
                        data = y,
                        burnin = 1000, 
                        pdf_func = pdf_func,
                        dist_type = "continuous",
                        loc = "mu",
                        vars_to_keep = c("eta", "mu", "sigma")) #vars_to_keep ignore numbers

# plot the mixture
plot(bmix)
```

![](external_comp_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# mode estimation
bayesmode = bayes_mode(bmix)

# plot mode inference
plot(bayesmode)
```

![](external_comp_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
# Summary of mode inference
summary(bayesmode)
```

    ## Posterior probability of multimodality is 0.058 
    ## 
    ## Snapshot of inference results on the number of modes:
    ##   tb_nb_modes (matrix, dim 3x2): 
    ##      number of modes posterior probability
    ## [1,]               1                 0.942
    ## [2,]               2                 0.046
    ## [3,]               3                 0.012
    ## 
    ## Snapshot of inference results on mode locations:
    ##   table_location (matrix, dim 75x2): 
    ##      mode location posterior probability
    ## [1,]           4.7                 0.001
    ## [2,]           4.8                 0.002
    ## [3,]           4.9                 0.010
    ## [4,]           5.0                 0.049
    ## [5,]           5.1                 0.000
    ## [6,]           5.2                 0.125
    ## ... (69 more rows)

### BNPmix

#### Estimation

``` r
set.seed(123)
library(BNPmix)
library(dplyr)

y = c(rnorm(100,0,1),
      rnorm(100,5,1))

## estimation
fit = PYdensity(y,
                mcmc = list(niter = 2000,
                            nburn = 1000,
                            print_message = FALSE),
                output = list(out_param = TRUE))
```

#### Transforming the output into a mcmc matrix with one column per variable

``` r
fit_mat = list()

for (i in 1:length(fit$p)) {
  k = length(fit$p[[i]][, 1])
  
  draw = c(fit$p[[i]][, 1],
           fit$mean[[i]][, 1],
           sqrt(fit$sigma2[[i]][, 1]),
           i)
  
  names(draw)[1:k] = paste0("eta", 1:k)
  names(draw)[(k+1):(2*k)] = paste0("mu", 1:k)
  names(draw)[(2*k+1):(3*k)] = paste0("sigma", 1:k)
  names(draw)[3*k + 1] = "draw"
  
  fit_mat[[i]] = draw
}

fit_mat = as.matrix(bind_rows(fit_mat))
```

#### Create a BayesMixture object

``` r
bmix = bayes_mixture(mcmc = fit_mat,
                        data = y,
                        burnin = 0, # the burnin has already been discarded
                        dist = "normal",
                        vars_to_keep = c("eta", "mu", "sigma"))

# plot the mixture
plot(bmix)
```

![](external_comp_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# mode estimation
bayesmode = bayes_mode(bmix)

# plot mode inference
plot(bayesmode)
```

![](external_comp_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# Summary of mode inference
summary(bayesmode)
```

    ## Posterior probability of multimodality is 1 
    ## 
    ## Snapshot of inference results on the number of modes:
    ##   tb_nb_modes (matrix, dim 2x2): 
    ##      number of modes posterior probability
    ## [1,]               2                 0.993
    ## [2,]               3                 0.007
    ## 
    ## Snapshot of inference results on mode locations:
    ##   table_location (matrix, dim 85x2): 
    ##      mode location posterior probability
    ## [1,] -3.000000e-01                 0.006
    ## [2,] -2.000000e-01                 0.000
    ## [3,] -1.000000e-01                 0.000
    ## [4,]  5.551115e-17                 0.000
    ## [5,]  1.000000e-01                 0.000
    ## [6,]  2.000000e-01                 0.192
    ## ... (79 more rows)

## Mode estimation in mixtures estimated with maximum likelihood

In the examples presented below external R package are used for
estimating mixtures with the EM algorithm while `BayesMultiMode` is used
for estimating and plotting modes.

### mixtools

``` r
set.seed(123)
library(mixtools)

y = c(rnorm(100,0,1),
      rnorm(100,3.5,1.5))

fit = normalmixEM(y)
```

    ## number of iterations= 172

``` r
pars = c(eta = fit$lambda, mu = fit$mu, sigma = fit$sigma)

mix = mixture(pars, dist = "normal", range = c(min(y), max(y))) # create a new object of class Mixture
modes = mix_mode(mix) # estimate modes and create an object of class Mode

summary(modes)
```

    ## Modes of a normal mixture with 2 components.
    ## - Number of modes found: 2
    ## - Mode estimation technique: fixed-point algorithm
    ## - Estimates of mode locations:
    ##   mode_estimates (numeric vector, dim 2): 
    ## [1] 3 0

``` r
# mode inference
plot(modes)
```

![](external_comp_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
summary(modes)
```

    ## Modes of a normal mixture with 2 components.
    ## - Number of modes found: 2
    ## - Mode estimation technique: fixed-point algorithm
    ## - Estimates of mode locations:
    ##   mode_estimates (numeric vector, dim 2): 
    ## [1] 3 0

### mclust

``` r
set.seed(123)
library(mclust)

y = c(rnorm(100,0,1),
      rnorm(100,3.5,1.5))

fit = Mclust(y)

pars = c(eta = fit$parameters$pro,
         mu = fit$parameters$mean,
         sigma = sqrt(fit$parameters$variance$sigmasq))

mix = mixture(pars, dist = "normal", range = c(min(y), max(y))) # creates a new object of class Mixture

summary(mix)
```

    ## 
    ##  Estimated mixture distribution.
    ## - Mixture type: continuous
    ## - Number of components: 2
    ## - Distribution family: normal
    ## - Number of distribution variables: 2
    ## - Names of variables: mu sigma
    ## - Parameter estimates:
    ##   pars (numeric vector, dim 6): 
    ##        eta1        eta2        mu.1        mu.2      sigma1      sigma2 
    ##  0.38684421  0.61315579 -0.03856849  2.82059287  0.83331314  1.75291009

``` r
# mode inference
modes = mix_mode(mix) # estimates modes and creates an object of class Mode

plot(modes)
```

![](external_comp_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
summary(modes)
```

    ## Modes of a normal mixture with 2 components.
    ## - Number of modes found: 2
    ## - Mode estimation technique: fixed-point algorithm
    ## - Estimates of mode locations:
    ##   mode_estimates (numeric vector, dim 2): 
    ## [1] 0 3
