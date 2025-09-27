# bayes_mode works with external MCMC output

    Code
      summary(bayesmix)
    Output
      Mixture estimated with a Bayesian MCMC method.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: NA
      - Number of distribution variables: 4
      - Names of variables: mu sigma xi nu
      
      Summary of MCMC output after burnin:
      # A tibble: 10 x 10
         variable  mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
         <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
       1 eta1       0.8    0.8    NA     0   0.8   0.8    NA       NA       NA
       2 eta2       0.2    0.2    NA     0   0.2   0.2    NA       NA       NA
       3 mu1        0.5    0.5    NA     0   0.5   0.5    NA       NA       NA
       4 mu2        6      6      NA     0   6     6      NA       NA       NA
       5 sigma1     1      1      NA     0   1     1      NA       NA       NA
       6 sigma2     2      2      NA     0   2     2      NA       NA       NA
       7 xi1        0      0      NA     0   0     0      NA       NA       NA
       8 xi2        0      0      NA     0   0     0      NA       NA       NA
       9 nu1        5      5      NA     0   5     5      NA       NA       NA
      10 nu2        5      5      NA     0   5     5      NA       NA       NA
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      Posterior probability of multimodality is 1 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 55x2): 
           mode location posterior probability
      [1,]           0.5                     1
      [2,]           0.6                     0
      [3,]           0.7                     0
      [4,]           0.8                     0
      [5,]           0.9                     0
      [6,]           1.0                     0
      ... (49 more rows)

# bayes_mode works with normal mixture

    Code
      summary(bayesmix)
    Output
      Mixture estimated with a Bayesian MCMC method.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: normal
      - Number of distribution variables: 2
      - Names of variables: mu sigma
      
      Summary of MCMC output after burnin:
      # A tibble: 6 x 10
        variable   mean median     sd    mad     q5    q95  rhat ess_bulk ess_tail
        <chr>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
      1 eta1      0.502  0.502 0.0366 0.0347  0.439  0.560 1.00      423.     361.
      2 eta2      0.498  0.498 0.0366 0.0347  0.440  0.561 1.00      423.     361.
      3 mu1       4.98   4.98  0.108  0.113   4.81   5.17  1.01      288.     288.
      4 mu2      -5.00  -5.00  0.0911 0.0914 -5.15  -4.85  0.998     353.     404.
      5 sigma1    0.998  0.998 0.0693 0.0730  0.887  1.12  1.00      325.     337.
      6 sigma2    0.907  0.904 0.0610 0.0595  0.814  1.02  1.00      303.     312.
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      Posterior probability of multimodality is 1 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 107x2): 
           mode location posterior probability
      [1,]          -5.2                 0.055
      [2,]          -5.1                 0.000
      [3,]          -5.0                 0.415
      [4,]          -4.9                 0.225
      [5,]          -4.8                 0.045
      [6,]          -4.7                 0.000
      ... (101 more rows)

# bayes_mode works with skew_normal mixture

    Code
      summary(bayesmix)
    Output
      Mixture estimated with a Bayesian MCMC method.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: skew_normal
      - Number of distribution variables: 3
      - Names of variables: xi omega alpha
      
      Summary of MCMC output after burnin:
      # A tibble: 8 x 10
        variable   mean median     sd    mad     q5    q95  rhat ess_bulk ess_tail
        <chr>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
      1 eta1      0.502  0.503 0.0350 0.0343  0.442  0.557 1.01     318.     363. 
      2 eta2      0.498  0.497 0.0350 0.0343  0.443  0.558 1.01     318.     363. 
      3 xi1       3.91   3.91  0.153  0.153   3.64   4.17  0.999     93.2    190. 
      4 xi2      -4.20  -4.17  0.245  0.233  -4.66  -3.87  1.02      42.2     80.2
      5 omega1    1.46   1.46  0.150  0.161   1.25   1.72  0.999    108.     121. 
      6 omega2    1.22   1.22  0.170  0.178   0.932  1.49  1.02      63.7    109. 
      7 alpha1    2.16   2.11  0.540  0.491   1.46   3.19  1.00      54.1     78.9
      8 alpha2   -1.39  -1.41  0.565  0.583  -2.34  -0.434 1.04      32.5     51.5
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      Posterior probability of multimodality is 1 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 103x2): 
           mode location posterior probability
      [1,]          -5.2                0.0025
      [2,]          -5.1                0.0000
      [3,]          -5.0                0.1325
      [4,]          -4.9                0.2750
      [5,]          -4.8                0.2675
      [6,]          -4.7                0.1900
      ... (97 more rows)

# bayes_mode works with shifted poisson mixture

    Code
      summary(bayesmix)
    Output
      Mixture estimated with a Bayesian MCMC method.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: shifted_poisson
      - Number of distribution variables: 2
      - Names of variables: kappa lambda
      
      Summary of MCMC output after burnin:
      # A tibble: 6 x 10
        variable  mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
        <chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
      1 eta1     0.435  0.435 0.0518 0.0548 0.352 0.523  1.20     3.81     63.9
      2 eta2     0.565  0.565 0.0518 0.0548 0.477 0.648  1.20     3.81     63.9
      3 kappa1   0      0     0      0      0     0     NA       NA        NA  
      4 kappa2   1.14   1     1.28   1.48   0     3      1.43     2.11     NA  
      5 lambda1  0.828  0.825 0.128  0.123  0.621 1.04   1.12     6.76    182. 
      6 lambda2  4.35   4.60  1.08   1.10   2.67  5.66   1.29     3.15     43.4
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      Posterior probability of multimodality is 0.9975 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 2x2): 
           number of modes posterior probability
      [1,]               1                0.0025
      [2,]               2                0.9975
      
      Inference results on mode locations:
        p_loc (matrix, dim 7x2): 
           mode location posterior probability
      [1,]             0                0.8525
      [2,]             1                0.1475
      [3,]             2                0.0000
      [4,]             3                0.0000
      [5,]             4                0.1325
      [6,]             5                0.8275
      ... (1 more rows)

# bayes_mode works with poisson mixture

    Code
      summary(bayesmix)
    Output
      Mixture estimated with a Bayesian MCMC method.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: poisson
      - Number of distribution variables: 1
      - Names of variables: lambda
      
      Summary of MCMC output after burnin:
      # A tibble: 4 x 10
        variable   mean median     sd    mad    q5    q95  rhat ess_bulk ess_tail
        <chr>     <dbl>  <dbl>  <dbl>  <dbl> <dbl>  <dbl> <dbl>    <dbl>    <dbl>
      1 eta1      0.503  0.503 0.0332 0.0319 0.447  0.556 0.998     359.     333.
      2 eta2      0.497  0.497 0.0332 0.0319 0.444  0.553 0.998     359.     333.
      3 lambda1   0.483  0.479 0.0721 0.0731 0.372  0.605 1.00      338.     271.
      4 lambda2  10.1   10.1   0.316  0.320  9.61  10.6   1.00      416.     452.
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      Posterior probability of multimodality is 1 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 12x2): 
           mode location posterior probability
      [1,]             0                     1
      [2,]             1                     0
      [3,]             2                     0
      [4,]             3                     0
      [5,]             4                     0
      [6,]             5                     0
      ... (6 more rows)

