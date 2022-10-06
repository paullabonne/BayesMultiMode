BayesMultiMode
================

This package provides functions for estimating the number of modes and
their associated locations in mixtures of discrete distributions using
Bayesian methods. The testing approach works in three stages. First, a
mixture distribution is estimated on the data using MCMC methods. The
number of mixture components is allowed to be unknown a priori and
estimated using a Spare Finite Mixture (SFM) algorithm. Second, a simple
post-processing is applied to the MCMC output where a given number of
iterations are used for burn-in and empty components are discarded.
Third, post-processed MCMC draws are used to compute the number of modes
as well as their locations and posterior probabilities.

Currently the package supports a flexible mixture of shifted poisson
distributions. The shifted poisson augments the poisson distribution
with a location parameter. More mixtures of discrete and continuous
distributions are in the pipeline.

#### Installing BayesMultiMode

``` r
# Installing from github :
# install.packages("devtools") # if devtools is not installed 
devtools::install_github("paullabonne/BayesMultiMode")
```

#### Loading the package

``` r
library(BayesMultiMode)
```

#### Generating data

``` r
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

#### Choosing either simulated or DNA data

``` r
# Select DNA data :
data("d4z4")
y = d4z4

# Or select simulated data :
# y = simulated_data
```

#### Setting parameters for SFM MCMC estimation

``` r
# Number of MCMC iterations 
M = 5000 

# Proportion of draws to discard as burnin
S = 0.5 

# Maximum number of mixture components 
Jmax = 6
```

#### Estimation with SFM MCMC

``` r
#Bayesian estimation
sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
```

    ## 10  % draws finished. Accept. prob of e0 = 9 percent
    ## 20  % draws finished. Accept. prob of e0 = 8 percent
    ## 30  % draws finished. Accept. prob of e0 = 8 percent
    ## 40  % draws finished. Accept. prob of e0 = 9 percent
    ## 50  % draws finished. Accept. prob of e0 = 9 percent
    ## 60  % draws finished. Accept. prob of e0 = 9 percent
    ## 70  % draws finished. Accept. prob of e0 = 9 percent
    ## 80  % draws finished. Accept. prob of e0 = 9 percent
    ## 90  % draws finished. Accept. prob of e0 = 9 percent
    ## 100  % draws finished. Accept. prob of e0 = 9 percent

``` r
#Plots of the estimation output
graphs = plots_mcmc(sfm_mcmc,S)
graphs[[3]]
```

<img src="README_files/figure-gfm/unnamed-chunk-6-1.png" width="70%" style="display: block; margin: auto;" />

``` r
graphs[[4]]
```

    ## Warning: Removed 41 rows containing missing values (position_stack).

<img src="README_files/figure-gfm/unnamed-chunk-6-2.png" width="70%" style="display: block; margin: auto;" />

#### Post-processing : burn-in and discarding empty components

``` r
post_sfmmcmc = post_sfm_mcmc(sfm_mcmc,S)
```

#### Mode inference

``` r
sfm_mcmc_modes = bayes_mode(post_sfmmcmc$theta_draws_slim,y)
sfm_mcmc_modes$graphs
```

<img src="README_files/figure-gfm/unnamed-chunk-8-1.png" width="70%" style="display: block; margin: auto;" />
