# mix_mode() returns expected results with dist = shifted_poisson and flat modes

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: shifted_poisson
      - Number of distribution variables: 2
      - Names of variables: lambda kappa

---

    Code
      summary(modes)
    Output
      
       Modes of a shifted_poisson mixture with 2 components.
      - Number of modes found: 3
      - Mode estimation technique: discrete algorithm

---

    Code
      modes$mode_estimates
    Output
      [1]  0  1 10

# mix_mode() function returns expected results with dist = poisson

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: poisson
      - Number of distribution variables: 1
      - Names of variables: lambda

---

    Code
      summary(modes)
    Output
      
       Modes of a poisson mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: discrete algorithm

---

    Code
      modes$mode_estimates
    Output
      [1] 0 9

# mix_mode() function returns expected results with arbitrary function

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: discrete
      - Number of components: 2
      - Distribution family: NA
      - Number of distribution variables: 2
      - Names of variables: mu size

---

    Code
      summary(modes)
    Output
      
       Modes of a discrete mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: discrete algorithm

---

    Code
      modes$mode_estimates
    Output
      [1]  0 18

# mix_mode() function returns expected results with dist = skew_normal

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: skew_normal
      - Number of distribution variables: 3
      - Names of variables: xi omega alpha

---

    Code
      summary(modes)
    Output
      
       Modes of a skew_normal mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm

---

    Code
      modes$mode_estimates
    Output
      [1] 0.002088749 5.999997076

# mix_mode() function returns expected results with an arbitrary function

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: NA
      - Number of distribution variables: 4
      - Names of variables: mu sigma xi nu

---

    Code
      summary(modes)
    Output
      
       Modes of a continuous mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm

---

    Code
      modes$mode_estimates
    Output
      [1] 0.00182144 5.88332478

# mix_mode() function returns expected results

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
      - Mixture type: continuous
      - Number of components: 2
      - Distribution family: normal
      - Number of distribution variables: 2
      - Names of variables: mu sigma

---

    Code
      summary(modes)
    Output
      
       Modes of a normal mixture with 2 components.
      - Number of modes found: 2
      - Mode estimation technique: fixed-point algorithm

---

    Code
      modes$mode_estimates
    Output
      [1] 0.006915293 4.999402022

