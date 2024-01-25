BayesMultiMode
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/paullabonne/BayesMultiMode/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paullabonne/BayesMultiMode/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/BayesMultiMode)](https://cran.r-project.org/package=BayesMultiMode)
<!-- badges: end -->

`BayesMultiMode` is an R package for detecting and exploring
multimodality using Bayesian techniques. The approach works in two
stages. First, a mixture distribution is fitted on the data using a
sparse finite mixture Markov chain Monte Carlo (SFM MCMC) algorithm. The
number of mixture components does not have to be known; the size of the
mixture is estimated endogenously through the SFM approach. Second, the
modes of the estimated mixture in each MCMC draw are retrieved using
algorithms specifically tailored for mode detection. These estimates are
then used to construct posterior probabilities for the number of modes,
their locations and uncertainties, providing a powerful tool for mode
inference. See Basturk et al. (2023) and Cross et al. (2023) for more
details.

### Installing BayesMultiMode from CRAN

``` r
install.packages("BayesMultiMode")
```

### Or installing the development version from GitHub

``` r
# install.packages("devtools") # if devtools is not installed 
devtools::install_github("paullabonne/BayesMultiMode")
```

### Loading BayesMultiMode

``` r
library(BayesMultiMode)
```

### Using BayesMultiMode for both MCMC estimation and mode inference

`BayesMultiMode` provides a very flexible and efficient MCMC estimation
approach : it handles mixtures with unknown number of components through
the sparse finite mixture approach of Malsiner-Walli,
Fruhwirth-Schnatter, and Grun (2016) and supports a comprehensive range
of mixture distributions, both continuous and discrete.

#### Estimation

``` r
set.seed(123)

# retrieve galaxy data
y = galaxy

# estimation
bayesmix = bayes_estimation(data = y,
                            K = 10,
                            dist = "normal",
                            nb_iter = 2000,
                            burnin = 1000,
                            printing = F)

plot(bayesmix, draws = 200)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="70%" style="display: block; margin: auto;" />

``` r
summary(bayesmix)
```

    ## 
    ##  Mixture estimated with a Bayesian MCMC method.
    ## - Mixture type: continuous
    ## - Number of components: 10
    ## - Distribution family: normal
    ## - Number of distribution variables: 2
    ## - Names of variables: mu sigma
    ## 
    ## Summary of MCMC output after burnin:

    ## # A tibble: 30 × 10
    ##    variable     mean   median      sd      mad        q5      q95  rhat ess_bulk
    ##    <chr>       <num>    <num>   <num>    <num>     <num>    <num> <num>    <num>
    ##  1 eta1     0.000214 2.50e-34 0.00224 3.71e-34 5.09e-272  6.32e-5 1.01     21.7 
    ##  2 eta2     0.00589  2.28e-29 0.0183  3.39e-29 1.73e-302  4.21e-2 1.04     12.4 
    ##  3 eta3     0.445    3.82e- 1 0.416   5.66e- 1 1.54e- 64  9.16e-1 2.12      1.31
    ##  4 eta4     0.0864   8.27e- 2 0.0313  3.05e- 2 4.31e-  2  1.44e-1 0.999  1009.  
    ##  5 eta5     0.000202 6.55e-34 0.00153 9.71e-34 7.96e-312  1.03e-4 1.02     18.6 
    ##  6 eta6     0.000300 6.43e-36 0.00245 9.53e-36 8.04e-287  2.38e-4 1.03     19.3 
    ##  7 eta7     0.436    5.08e- 1 0.424   6.18e- 1 1.68e-267  9.13e-1 2.15      1.30
    ##  8 eta8     0.00932  1.25e-16 0.0203  1.85e-16 2.22e-109  5.64e-2 1.30      2.68
    ##  9 eta9     0.00677  2.63e-25 0.0172  3.89e-25 1.07e-242  4.66e-2 1.05     15.6 
    ## 10 eta10    0.00977  2.83e-19 0.0189  4.19e-19 1.57e-245  5.39e-2 1.30      2.76
    ## # ℹ 20 more rows
    ## # ℹ 1 more variable: ess_tail <num>

#### Mode inference

``` r
# mode estimation
bayesmode = bayes_mode(bayesmix)

plot(bayesmode)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="70%" style="display: block; margin: auto;" />

``` r
summary(bayesmode)
```

    ## The posterior probability of the data being multimodal is 0.999
    ## 
    ##  Number of estimated modes and their posterior probabilities:

    ##      Number of modes Posterior probabilty
    ## [1,]               1                0.001
    ## [2,]               2                0.067
    ## [3,]               3                0.902
    ## [4,]               4                0.030

### Using BayesMultiMode for mode inference with external MCMC output

`BayesMultiMode` also works on MCMC output generated using external
software. The function `new_BayesMixture()` creates an object of class
`BayesMixture` which can then be used as input in the mode inference
function `bayes_mode()`. Here is an example using cyclone intensity data
(Knapp et al. 2018) and the `BNPmix` package for estimation.

``` r
library(BNPmix)
library(dplyr)

y = cyclone %>%
  filter(BASIN == "SI",
         SEASON > "1981") %>%
  select(max_wind) %>%
  unlist()

## estimation
PY_result = PYdensity(y,
                      mcmc = list(niter = 2000,
                                  nburn = 1000,
                                  print_message = FALSE),
                      output = list(out_param = TRUE))
```

#### Transforming the output into a mcmc matrix with one column per variable

``` r
library(dplyr)

mcmc_py = list()

for (i in 1:length(PY_result$p)) {
  k = length(PY_result$p[[i]][, 1])
  
  draw = c(PY_result$p[[i]][, 1],
           PY_result$mean[[i]][, 1],
           sqrt(PY_result$sigma2[[i]][, 1]),
           i)
  
  names(draw)[1:k] = paste0("eta", 1:k)
  names(draw)[(k+1):(2*k)] = paste0("mu", 1:k)
  names(draw)[(2*k+1):(3*k)] = paste0("sigma", 1:k)
  names(draw)[3*k + 1] = "draw"
  
  mcmc_py[[i]] = draw
}

mcmc_py = as.matrix(bind_rows(mcmc_py))
```

#### Creating an object of class `BayesMixture`

``` r
py_BayesMix = new_BayesMixture(mcmc = mcmc_py,
                               data = y,
                               K = (ncol(mcmc_py)-1)/3,
                               burnin = 0, # the burnin has already been discarded
                               dist = "normal",
                               dist_type = "continuous",
                               vars_to_keep = c("eta", "mu", "sigma"))

summary(py_BayesMix)
```

    ## 
    ##  Mixture estimated with a Bayesian MCMC method.
    ## - Mixture type: continuous
    ## - Number of components: 15
    ## - Distribution family: normal
    ## - Number of distribution variables: 2
    ## - Names of variables: mu sigma
    ## 
    ## Summary of MCMC output after burnin:

    ## # A tibble: 45 × 10
    ##    variable    mean  median      sd     mad        q5    q95  rhat ess_bulk
    ##    <chr>      <num>   <num>   <num>   <num>     <num>  <num> <num>    <num>
    ##  1 eta1      0.292   0.299   0.0483  0.0405  0.197     0.357  1.05     20.3
    ##  2 eta2      0.285   0.286   0.0546  0.0612  0.194     0.370  1.03     27.1
    ##  3 eta3      0.368   0.367   0.0455  0.0488  0.298     0.443  1.09     29.7
    ##  4 eta4      0.0412  0.0212  0.0488  0.0297  0.000423  0.138  1.03     17.1
    ##  5 eta5     NA      NA      NA      NA      NA        NA     NA        NA  
    ##  6 eta6     NA      NA      NA      NA      NA        NA     NA        NA  
    ##  7 eta7     NA      NA      NA      NA      NA        NA     NA        NA  
    ##  8 eta8     NA      NA      NA      NA      NA        NA     NA        NA  
    ##  9 eta9     NA      NA      NA      NA      NA        NA     NA        NA  
    ## 10 eta10    NA      NA      NA      NA      NA        NA     NA        NA  
    ## # ℹ 35 more rows
    ## # ℹ 1 more variable: ess_tail <num>

#### Plotting the mixture

``` r
plot(py_BayesMix)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="70%" style="display: block; margin: auto;" />

#### Mode inference

``` r
# mode estimation
bayesmode = bayes_mode(py_BayesMix)

# plot
plot(bayesmode)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="70%" style="display: block; margin: auto;" />

``` r
# Summary 
summary(bayesmode)
```

    ## The posterior probability of the data being multimodal is 1
    ## 
    ##  Number of estimated modes and their posterior probabilities:

    ##      Number of modes Posterior probabilty
    ## [1,]               2                0.897
    ## [2,]               3                0.102
    ## [3,]               4                0.001

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-basturk_2023" class="csl-entry">

Basturk, Nalan, Jamie L. Cross, Peter de Knijff, Lennart Hoogerheide,
Paul Labonne, and Herman K. van Dijk. 2023. “BayesMultiMode: Bayesian
Mode Inference in r.” *Tinbergen Institute Discussion Paper TI
2023-041/III*.

</div>

<div id="ref-cross_2023" class="csl-entry">

Cross, Jamie L., Lennart Hoogerheide, Paul Labonne, and Herman K. van
Dijk. 2023. “Credible Mode Determination in Multimodal Economic and
Financial Data Distributions.” *Tinbergen Institute Discussion Paper TI
2023-038/III*.

</div>

<div id="ref-knapp_international_2018" class="csl-entry">

Knapp, Kenneth R., Howard J. Diamond, Kossin J. P., Michael C. Kruk, and
C. J. Schreck. 2018. “International Best Track Archive for Climate
Stewardship (IBTrACS) Project, Version 4.” *NOAA National Centers for
Environmental Information*. <https://doi.org/10.1175/2009BAMS2755.1>.

</div>

<div id="ref-malsiner-walli_model-based_2016" class="csl-entry">

Malsiner-Walli, Gertraud, Sylvia Fruhwirth-Schnatter, and Bettina Grun.
2016. “Model-Based Clustering Based on Sparse Finite Gaussian Mixtures.”
*Statistics and Computing* 26 (1): 303–24.
<https://doi.org/10.1007/s11222-014-9500-2>.

</div>

</div>
