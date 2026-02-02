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
      The posterior probability of multimodality is 1 
      
       The most likely number of modes is 2 
      
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
      1 eta1      0.499  0.500 0.0347 0.0379  0.444  0.552 1.000     369.     373.
      2 eta2      0.501  0.500 0.0347 0.0379  0.448  0.556 1.000     369.     373.
      3 mu1      -4.99  -4.99  0.0836 0.0900 -5.13  -4.86  1.00      384.     371.
      4 mu2       4.98   4.98  0.0955 0.0908  4.82   5.14  1.00      380.     389.
      5 sigma1    0.915  0.914 0.0588 0.0572  0.819  1.01  0.998     444.     303.
      6 sigma2    0.993  0.995 0.0685 0.0688  0.883  1.10  1.00      404.     365.
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of multimodality is 1 
      
       The most likely number of modes is 2 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 106x2): 
           mode location posterior probability
      [1,]          -5.2                0.0375
      [2,]          -5.1                0.0000
      [3,]          -5.0                0.4250
      [4,]          -4.9                0.2925
      [5,]          -4.8                0.0325
      [6,]          -4.7                0.0025
      ... (100 more rows)

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
        variable   mean median     sd    mad      q5    q95  rhat ess_bulk ess_tail
        <chr>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>    <dbl>    <dbl>
      1 eta1      0.497  0.497 0.0348 0.0347  0.440   0.553 1.000    344.     365. 
      2 eta2      0.503  0.503 0.0348 0.0347  0.447   0.560 1.000    344.     365. 
      3 xi1      -4.19  -4.16  0.228  0.213  -4.57   -3.85  1.03      37.4    125. 
      4 xi2       4.44   4.40  0.309  0.301   3.99    4.96  1.02      40.8    103. 
      5 omega1    1.23   1.22  0.161  0.155   0.972   1.49  1.03      42.5    119. 
      6 omega2    1.24   1.23  0.148  0.142   1.02    1.53  1.00      66.2     80.8
      7 alpha1   -1.34  -1.27  0.526  0.488  -2.25   -0.553 1.06      25.6     45.8
      8 alpha2    0.634  0.634 0.421  0.406  -0.0398  1.29  1.02      32.2     90.7
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of multimodality is 1 
      
       The most likely number of modes is 2 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 1x2): 
           number of modes posterior probability
      [1,]               2                     1
      
      Inference results on mode locations:
        p_loc (matrix, dim 106x2): 
           mode location posterior probability
      [1,]          -5.2                0.0025
      [2,]          -5.1                0.0000
      [3,]          -5.0                0.0900
      [4,]          -4.9                0.2750
      [5,]          -4.8                0.2975
      [6,]          -4.7                0.1975
      ... (100 more rows)

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
      1 eta1     0.417  0.417 0.0492 0.0498 0.336 0.498  1.02    44.4     102. 
      2 eta2     0.583  0.583 0.0492 0.0498 0.502 0.664  1.02    44.4     102. 
      3 kappa1   0      0     0      0      0     0     NA       NA        NA  
      4 kappa2   0.650  0     0.954  0      0     3      1.08     8.49     NA  
      5 lambda1  0.791  0.794 0.123  0.124  0.613 0.994  1.01    54.9     230. 
      6 lambda2  4.71   4.98  0.826  0.686  2.84  5.68   1.05    14.9      11.9
      this table can be reproduced with: summarise_draws(bayesmix$mcmc)
      
      Note that label-switching might occur in the MCMC draws becayse BayesMultiMode does not carry out post-processing. 
      While label-switching does not affect mode inference it can affect diagnostic checks.
    Message
      

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of multimodality is 0.9975 
      
       The most likely number of modes is 2 
      
      Inference results on the number of modes:
        p_nb_modes (matrix, dim 2x2): 
           number of modes posterior probability
      [1,]               1                0.0025
      [2,]               2                0.9975
      
      Inference results on mode locations:
        p_loc (matrix, dim 7x2): 
           mode location posterior probability
      [1,]             0                0.9025
      [2,]             1                0.0975
      [3,]             2                0.0000
      [4,]             3                0.0000
      [5,]             4                0.2025
      [6,]             5                0.7750
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
      The posterior probability of multimodality is 1 
      
       The most likely number of modes is 2 
      
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

