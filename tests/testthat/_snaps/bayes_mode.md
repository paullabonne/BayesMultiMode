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
        p_loc (matrix, dim 104x2): 
           mode location posterior probability
      [1,]          -5.2                0.0050
      [2,]          -5.1                0.0000
      [3,]          -5.0                0.2800
      [4,]          -4.9                0.3600
      [5,]          -4.8                0.2075
      [6,]          -4.7                0.0675
      ... (98 more rows)

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

