# mix_mode() returns expected results with dist = shifted_poisson and flat modes

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
       Number of components: 2
       Number of component parameters: 2
       Mixture family: shifted_poisson
       Mixture type: discrete

---

    Code
      summary(modes)
    Output
      
       Modes of a shifted_poisson mixture with 2 components.
       Number of modes found: 3
       Mode estimation technique: discrete algorithm

---

    Code
      modes$mode_estimates
    Output
      [1]  0  1 10

# discrete_MF function returns expected results with dist = poisson

    Code
      summary(mix)
    Output
      
       Estimated mixture distribution.
       Number of components: 2
       Number of component parameters: 1
       Mixture family: poisson
       Mixture type: discrete

---

    Code
      summary(modes)
    Output
      
       Modes of a poisson mixture with 2 components.
       Number of modes found: 2
       Mode estimation technique: discrete algorithm

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
       Number of components: 2
       Number of component parameters: 2
       Mixture family: NA
       Mixture type: discrete

---

    Code
      summary(modes)
    Output
      
       Modes of a discrete mixture with 2 components.
       Number of modes found: 2
       Mode estimation technique: discrete algorithm

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
       Number of components: 2
       Number of component parameters: 3
       Mixture family: skew_normal
       Mixture type: continuous

---

    Code
      summary(modes)
    Output
      
       Modes of a skew_normal mixture with 2 components.
       Number of modes found: 2
       Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm

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
       Number of components: 2
       Number of component parameters: 4
       Mixture family: NA
       Mixture type: continuous

---

    Code
      summary(modes)
    Output
      
       Modes of a continuous mixture with 2 components.
       Number of modes found: 2
       Mode estimation technique: Modal Expectation-Maximization (MEM) algorithm

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
       Number of components: 2
       Number of component parameters: 2
       Mixture family: normal
       Mixture type: continuous

---

    Code
      summary(modes)
    Output
      
       Modes of a normal mixture with 2 components.
       Number of modes found: 2
       Mode estimation technique: fixed-point algorithm

---

    Code
      modes$mode_estimates
    Output
      [1] 0.02829676 4.99985083

