# mix_mode() returns expected results with dist = shifted_poisson and flat modes

    Code
      modes$mode_estimates
    Output
      [1]  0  1 10

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: shifted_poisson
      - Number of distribution variables: 2
      - Names of variables: lambda kappa
      - Parameter estimates:
        pars (numeric vector, dim 6): 
         eta1    eta2 lambda1 lambda2  kappa1  kappa2 
          0.5     0.5     0.1     1.0    10.0     0.0 

---

    Code
      summary(modes)
    Output
      Modes of a shifted_poisson mixture with 2 components.
      - Number of modes found: 3
      - Mode estimation technique: discrete algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 3): 
      [1]  0  1 10

# mix_mode() function returns expected results with dist = poisson

    Code
      modes$mode_estimates
    Output
      [1] 0 9

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: poisson
      - Number of distribution variables: 1
      - Names of variables: lambda
      - Parameter estimates:
        pars (numeric vector, dim 4): 
         eta1    eta2 lambda1 lambda2 
          0.5     0.5     0.1    10.0 

---

    Code
      summary(modes)
    Output
      Modes of a poisson mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: discrete algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 2): 
      [1] 0 9

# mix_mode() function returns expected results with arbitrary function

    Code
      modes$mode_estimates
    Output
      [1]  0 18

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: NA
      - Number of distribution variables: 2
      - Names of variables: mu size
      - Parameter estimates:
        pars (numeric vector, dim 6): 
       eta1  eta2   mu1   mu2 size1 size2 
        0.5   0.5  20.0   5.0  20.0   0.5 

---

    Code
      summary(modes)
    Output
      Modes of a discrete mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: discrete algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 2): 
      [1]  0 18

# mix_mode() function returns expected results with dist = skew_normal

    Code
      modes$mode_estimates
    Output
      [1] 0.002088749 5.999997076

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: skew_normal
      - Number of distribution variables: 3
      - Names of variables: xi omega alpha
      - Parameter estimates:
        pars (numeric vector, dim 8): 
        eta1   eta2    xi1    xi2 omega1 omega2 
         0.8    0.2    0.0    6.0    1.0    2.0 
      ... (2 more elements)

---

    Code
      summary(modes)
    Output
      Modes of a skew_normal mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 2): 
      [1] 0 6

# mix_mode() function returns expected results with an arbitrary function

    Code
      modes$mode_estimates
    Output
      [1] 0.00182144 5.88332478

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: NA
      - Number of distribution variables: 4
      - Names of variables: mu sigma xi nu
      - Parameter estimates:
        pars (numeric vector, dim 10): 
        eta1   eta2    mu1    mu2 sigma1 sigma2 
         0.8    0.2    0.0    6.0    1.0    2.0 
      ... (4 more elements)

---

    Code
      summary(modes)
    Output
      Modes of a continuous mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 2): 
      [1] 0 6

# mix_mode() function returns expected results

    Code
      modes$mode_estimates
    Output
      [1] 0.006915293 4.999402022

---

    Code
      summary(mix)
    Output
      Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: normal
      - Number of distribution variables: 2
      - Names of variables: mu sigma
      - Parameter estimates:
        pars (numeric vector, dim 6): 
        eta1   eta2    mu1    mu2 sigma1 sigma2 
         0.8    0.2    0.0    5.0    1.0    2.0 

---

    Code
      summary(modes)
    Output
      Modes of a normal mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: fixed-point algorithm
      - Estimates of mode locations:
        mode_estimates (numeric vector, dim 2): 
      [1] 0 5

