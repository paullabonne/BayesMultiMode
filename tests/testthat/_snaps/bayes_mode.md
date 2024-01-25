# bayes_mode works with external MCMC output

    Code
      summary(bayesmix)
    Output
      
       Mixture estimated with a Bayesian MCMC method.
       Number of components: 2
       Number of component parameters: 4
       Mixture family: NA
       Mixture type: continuous

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of the data being multimodal is 1
      
       Number of estimated modes and their posterior probabilities:
           Number of modes Posterior probabilty 
                         2                    1 

# bayes_mode works with normal mixture

    Code
      summary(bayesmix)
    Output
      
       Mixture estimated with a Bayesian MCMC method.
       Number of components: 2
       Number of component parameters: 2
       Mixture family: normal
       Mixture type: continuous

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of the data being multimodal is 1
      
       Number of estimated modes and their posterior probabilities:
           Number of modes Posterior probabilty 
                         2                    1 

# bayes_mode works with skew_normal mixture

    Code
      summary(bayesmix)
    Output
      
       Mixture estimated with a Bayesian MCMC method.
       Number of components: 2
       Number of component parameters: 3
       Mixture family: skew_normal
       Mixture type: continuous

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of the data being multimodal is 1
      
       Number of estimated modes and their posterior probabilities:
           Number of modes Posterior probabilty 
                         2                    1 

# bayes_mode works with shifted poisson mixture

    Code
      summary(bayesmix)
    Output
      
       Mixture estimated with a Bayesian MCMC method.
       Number of components: 2
       Number of component parameters: 2
       Mixture family: shifted_poisson
       Mixture type: discrete

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of the data being multimodal is 0.9975
      
       Number of estimated modes and their posterior probabilities:
           Number of modes Posterior probabilty
      [1,]               1               0.0025
      [2,]               2               0.9975

# bayes_mode works with poisson mixture

    Code
      summary(bayesmix)
    Output
      
       Mixture estimated with a Bayesian MCMC method.
       Number of components: 2
       Number of component parameters: 1
       Mixture family: poisson
       Mixture type: discrete

---

    Code
      summary(bayesmode)
    Output
      The posterior probability of the data being multimodal is 1
      
       Number of estimated modes and their posterior probabilities:
           Number of modes Posterior probabilty 
                         2                    1 

