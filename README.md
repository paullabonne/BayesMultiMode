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
modes of the estimated mixture at each MCMC are retrieved using
algorithms specifically tailored for mode detection. These estimates are
then used to construct posterior probabilities for the number of modes,
their locations and uncertainties, providing a powerful tool for mode
inference.

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
approach : it handles mixtures with unknown number of components and
supports a comprehensive range of mixture distributions, both continuous
and discrete.

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
                            burnin = 1000)
```

    ## 10  % draws finished
    ## 20  % draws finished
    ## 30  % draws finished
    ## 40  % draws finished
    ## 50  % draws finished
    ## 60  % draws finished
    ## 70  % draws finished
    ## 80  % draws finished
    ## 90  % draws finished
    ## 100  % draws finished

``` r
# plot estimated mixture
plot(bayesmix, max_size = 200)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="70%" style="display: block; margin: auto;" />

#### Mode inference

``` r
# mode estimation
bayesmode = bayes_mode(bayesmix)

# plot 
plot(bayesmode, max_size = 200)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="70%" style="display: block; margin: auto;" />

``` r
# Summary 
summary(bayesmode)
```

    ## The posterior probability of the data being multimodal is 0.993 .
    ## 
    ## The number of estimated modes and their posterior probabilities is:

    ##      Number of modes Posterior probabilty
    ## [1,]               1                0.007
    ## [2,]               2                0.133
    ## [3,]               3                0.840
    ## [4,]               4                0.020

### Using BayesMultiMode for mode inference with external MCMC output

`BayesMultiMode` also works on MCMC output generated using external
software. The function `new_BayesMixture()` creates an object of class
`BayesMixture` which can then be used as input in the mode inference
function `bayes_mode()`. Here is an example where the `BNPmix` package
is used for MCMC estimation.

``` r
library(BNPmix)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
y = cyclone %>%
  filter(BASIN == "SI",
         SEASON > "1981") %>%
  select(max_wind) %>%
  unlist()

## estimation
PY_result = PYdensity(y,
                      mcmc = list(niter = 2000, nburn = 1000),
                      output = list(out_param = TRUE))
```

    ## Completed:   200/2000 - in 0.044338 sec
    ## Completed:   400/2000 - in 0.091791 sec
    ## Completed:   600/2000 - in 0.148233 sec
    ## Completed:   800/2000 - in 0.200009 sec
    ## Completed:   1000/2000 - in 0.248294 sec
    ## Completed:   1200/2000 - in 0.297524 sec
    ## Completed:   1400/2000 - in 0.349355 sec
    ## Completed:   1600/2000 - in 0.402126 sec
    ## Completed:   1800/2000 - in 0.457695 sec
    ## Completed:   2000/2000 - in 0.512106 sec
    ## 
    ## Estimation done in 0.512113 seconds

``` r
plot(PY_result)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="70%" style="display: block; margin: auto;" />

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

mcmc_py = bind_rows(mcmc_py)
```

#### Creating an object of class `BayesMixture`

``` r
pars_names = c(eta = "eta",
               mu = "mu",
               sigma = "sigma")

py_BayesMix = new_BayesMixture(mcmc = mcmc_py,
                               data = y,
                               K = (ncol(mcmc_py)-1)/3,
                               burnin = 0, # the burnin has already been discarded
                               dist = "normal",
                               pars_names = pars_names,
                               dist_type = "continuous")
```

#### Mode inference

``` r
# mode estimation
bayesmode = bayes_mode(py_BayesMix)

# plot ¨
plot(bayesmode, max_size = 200)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="70%" style="display: block; margin: auto;" />

``` r
# Summary 
summary(bayesmode)
```

    ## The posterior probability of the data being multimodal is 1 .
    ## 
    ## The number of estimated modes and their posterior probabilities is:

    ##      Number of modes Posterior probabilty
    ## [1,]               2                0.897
    ## [2,]               3                0.103

### References

(Malsiner-Walli, Fruhwirth-Schnatter, and Grun 2016) (Schaap et al.
2013)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-malsiner-walli_model-based_2016" class="csl-entry">

Malsiner-Walli, Gertraud, Sylvia Fruhwirth-Schnatter, and Bettina Grun.
2016. “Model-Based Clustering Based on Sparse Finite Gaussian Mixtures.”
*Statistics and Computing* 26 (1): 303–24.
<https://doi.org/10.1007/s11222-014-9500-2>.

</div>

<div id="ref-schaap_genome-wide_2013" class="csl-entry">

Schaap, Mireille, Richard JLF Lemmers, Roel Maassen, Patrick J. van der
Vliet, Lennart F. Hoogerheide, Herman K. van Dijk, Nalan Basturk, Peter
de Knijff, and Silvère M. van der Maarel. 2013. “Genome-Wide Analysis of
Macrosatellite Repeat Copy Number Variation in Worldwide Populations:
Evidence for Differences and Commonalities in Size Distributions and
Size Restrictions.” *BMC Genomics* 14 (1): 143.
<https://doi.org/10.1186/1471-2164-14-143>.

</div>

</div>
