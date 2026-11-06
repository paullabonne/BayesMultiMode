# Simulation designs: cases D1 to D4

These are the designs for the first parameter recovery check of the corrected shifted negative binomial mixture sampler.

## True parameters

r̃ = 1/(1 + r) and p = r/(r + μ).

| Case | n | K used | Component | Weight π | Shift κ | Mean μ | Shape r | r̃ | p | True modes |
|---|---|---|---|---|---|---|---|---|---|---|
| **D1** | 2,000 | 1 | 1 | 1.0 | 5 | 10 | 3 | 0.250 | 0.231 | 11 |
| **D2** | 10,000 | 2 | 1 | 0.5 | 0 | 15 | 10 | 0.091 | 0.400 | 13, 50 |
| | | | 2 | 0.5 | 50 | 1 | 2 | 0.333 | 0.667 | |
| **D3** | 2,000 | 2 | 1 | 0.6 | 0 | 2 | 5 | 0.167 | 0.714 | 1, 13 |
| | | | 2 | 0.4 | 10 | 5 | 3 | 0.250 | 0.375 | |
| **D4** | 2,000 | 6 | 1 | 0.6 | 0 | 2 | 5 | 0.167 | 0.714 | 1, 13 |
| | | | 2 | 0.4 | 10 | 5 | 3 | 0.250 | 0.375 | |
