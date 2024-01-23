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

---

    Code
      sum(bayesmode$modes, na.rm = T)
    Output
      [1] 6.426819

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

---

    Code
      sum(bayesmode$modes, na.rm = T)
    Output
      [1] -6.701141

---

    Code
      sum(bayesmix$loglik)
    Output
      [1] 14921.57

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

---

    Code
      sum(bayesmode$modes, na.rm = T)
    Output
      [1] -91.98267

---

    Code
      sum(bayesmix$loglik)
    Output
      [1] 14952.2

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

---

    Code
      sum(bayesmode$modes, na.rm = T)
    Output
      [1] 2016

---

    Code
      sum(bayesmix$loglik)
    Output
      [1] 12556.12

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

---

    Code
      sum(bayesmode$modes, na.rm = T)
    Output
      [1] 3855

---

    Code
      sum(bayesmix$loglik)
    Output
      [1] 14159.67

