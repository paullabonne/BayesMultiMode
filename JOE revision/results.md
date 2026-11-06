# D1-D4: single-chain simulation results

One simulated dataset and one MCMC chain per design. D3 and D4 use the same observations.

## Method

Each chain uses 10,000 iterations: 5,000 burn-in and 5,000 retained draws, without thinning. Fits use the current local shifted negative binomial sampler and its default priors. No chain was selected from multiple runs.

D1 calls the underlying sampler directly because the public mixture-object constructor rejects `K = 1`; D2-D4 use `bayes_fit`. All four use the same sampling algorithm and default prior settings (`a0 = 1`, `A0 = 200`, initial MH scale `sig = 0.1`).

The estimated distribution is the mean of the mixture PMF over all retained draws, including all six fitted components in D4. Shading gives 95% pointwise posterior intervals for the PMF, not predictive intervals for individual observations. Lines connect integer counts for readability.

Parameter estimates are posterior means, with equal-tailed 95% credible intervals. Shape and success probability are transformed in each draw before averaging: `r = 1/rtilde - 1` and `p = r/(r + mu)`. Here `mu` is the unshifted negative binomial mean; the component mean is `kappa + mu`.

For component comparisons, retain the largest-weight components (one in D1, two otherwise) and order them by `kappa + mu` in each draw. D4 extras are ranked by decreasing weight in each draw. This is a reporting convention, not evidence that exactly two components were recovered. Extras have true weight zero and undefined true shape/location parameters.

One chain and one dataset per design do not establish convergence or repeated-sampling recovery. No between-chain R-hat is available. Single-chain bulk ESS values below describe only the sampled trajectory.

## Distribution comparison

TV pools the probability above the evaluation grid into one tail bin; it is a lower bound on full-support total variation. CDF error is the maximum on the grid. Modes below belong to the true and posterior-mean PMFs, not the posterior distribution of mode counts.

| Case | n | K | True modes | Estimated modes | TV (tail pooled) | Max CDF error | Mean ESS | Log-likelihood ESS | Runtime (s) |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|
| D1 | 2000 | 1 | 11 | 11 | 0.0024 | 0.0018 | 3288 | 807 | 75.3 |
| D2 | 10000 | 2 | 13, 50 | 13, 50 | 0.0130 | 0.0108 | 3982 | 554 | 876.5 |
| D3 | 2000 | 2 | 1, 13 | 1, 13 | 0.0208 | 0.0120 | 3619 | 2 | 89.0 |
| D4 | 2000 | 6 | 1, 13 | 1, 13 | 0.0174 | 0.0106 | 4328 | 1045 | 103.1 |

D3 has the largest distribution error among these runs (TV 0.0208). The parameter tables below also show why a close distribution fit does not establish parameter recovery.

## D1

![D1 true and estimated probability mass functions](results-D1.png)

Data seed: `202601`. Chain seed: `31001`. PMF evaluation grid: 0-200. Estimated mass beyond the grid: 0.

Shift update: Joint shift-mean update (kappa + mu held fixed during the shift step).

| Component | Parameter | True | Estimated mean | 95% credible interval |
|---|---|---:|---:|---|
| 1 | Weight pi | 1.0000 | 1.0000 | [1.0000, 1.0000] |
| 1 | Shift kappa | 5.0000 | 5.0000 | [5.0000, 5.0000] |
| 1 | NB mean mu | 10.0000 | 9.9968 | [9.7077, 10.2870] |
| 1 | Shape r | 3.0000 | 2.9684 | [2.7522, 3.2083] |
| 1 | rtilde | 0.2500 | 0.2522 | [0.2376, 0.2665] |
| 1 | Success probability p | 0.2308 | 0.2289 | [0.2144, 0.2440] |
| 1 | Component mean kappa + mu | 15.0000 | 14.9968 | [14.7077, 15.2870] |

## D2

![D2 true and estimated probability mass functions](results-D2.png)

Data seed: `202602`. Chain seed: `31002`. PMF evaluation grid: 0-200. Estimated mass beyond the grid: 0.

Shift update: Joint shift-mean update (kappa + mu held fixed during the shift step).

| Component | Parameter | True | Estimated mean | 95% credible interval |
|---|---|---:|---:|---|
| 1 | Weight pi | 0.5000 | 0.5072 | [0.4976, 0.5168] |
| 1 | Shift kappa | 0.0000 | 0.0000 | [0.0000, 0.0000] |
| 1 | NB mean mu | 15.0000 | 14.8076 | [14.6370, 14.9765] |
| 1 | Shape r | 10.0000 | 9.4505 | [8.8624, 10.0661] |
| 1 | rtilde | 0.0909 | 0.0958 | [0.0904, 0.1014] |
| 1 | Success probability p | 0.4000 | 0.3895 | [0.3743, 0.4048] |
| 1 | Component mean kappa + mu | 15.0000 | 14.8076 | [14.6370, 14.9765] |
| 2 | Weight pi | 0.5000 | 0.4928 | [0.4832, 0.5024] |
| 2 | Shift kappa | 50.0000 | 50.0000 | [50.0000, 50.0000] |
| 2 | NB mean mu | 1.0000 | 1.0266 | [0.9925, 1.0623] |
| 2 | Shape r | 2.0000 | 1.9384 | [1.6916, 2.2202] |
| 2 | rtilde | 0.3333 | 0.3411 | [0.3105, 0.3715] |
| 2 | Success probability p | 0.6667 | 0.6530 | [0.6211, 0.6846] |
| 2 | Component mean kappa + mu | 51.0000 | 51.0266 | [50.9925, 51.0623] |

## D3

![D3 true and estimated probability mass functions](results-D3.png)

Data seed: `202603`. Chain seed: `31003`. PMF evaluation grid: 0-200. Estimated mass beyond the grid: 0.

Shift update: Joint shift-mean update (kappa + mu held fixed during the shift step).

**Mixing warning:** log-likelihood bulk ESS is only 2.4 out of 5000 retained draws. Treat the reported means and intervals as descriptions of this run, not reliable converged posterior summaries.

| Component | Parameter | True | Estimated mean | 95% credible interval |
|---|---|---:|---:|---|
| 1 | Weight pi | 0.6000 | 0.6103 | [0.5886, 0.6316] |
| 1 | Shift kappa | 0.0000 | 0.0000 | [0.0000, 0.0000] |
| 1 | NB mean mu | 2.0000 | 1.9993 | [1.9030, 2.1050] |
| 1 | Shape r | 5.0000 | 4.5814 | [3.4403, 6.2593] |
| 1 | rtilde | 0.1667 | 0.1820 | [0.1378, 0.2252] |
| 1 | Success probability p | 0.7143 | 0.6927 | [0.6306, 0.7625] |
| 1 | Component mean kappa + mu | 2.0000 | 1.9993 | [1.9030, 2.1050] |
| 2 | Weight pi | 0.4000 | 0.3897 | [0.3684, 0.4114] |
| 2 | Shift kappa | 10.0000 | 9.6794 | [9.0000, 10.0000] |
| 2 | NB mean mu | 5.0000 | 5.5380 | [4.9887, 6.3924] |
| 2 | Shape r | 3.0000 | 4.1085 | [2.8667, 6.6234] |
| 2 | rtilde | 0.2500 | 0.2048 | [0.1312, 0.2586] |
| 2 | Success probability p | 0.3750 | 0.4186 | [0.3525, 0.5150] |
| 2 | Component mean kappa + mu | 15.0000 | 15.2174 | [14.9660, 15.4912] |

## D4

![D4 true and estimated probability mass functions](results-D4.png)

Data seed: `202603`. Chain seed: `31004`. PMF evaluation grid: 0-200. Estimated mass beyond the grid: 1.4e-09.

Shift update: Joint shift-mean update (kappa + mu held fixed during the shift step).

| Component | Parameter | True | Estimated mean | 95% credible interval |
|---|---|---:|---:|---|
| 1 | Weight pi | 0.6000 | 0.6107 | [0.5894, 0.6316] |
| 1 | Shift kappa | 0.0000 | 0.0000 | [0.0000, 0.0000] |
| 1 | NB mean mu | 2.0000 | 2.0104 | [1.9073, 2.1128] |
| 1 | Shape r | 5.0000 | 4.5160 | [3.3782, 6.0163] |
| 1 | rtilde | 0.1667 | 0.1840 | [0.1425, 0.2284] |
| 1 | Success probability p | 0.7143 | 0.6887 | [0.6246, 0.7507] |
| 1 | Component mean kappa + mu | 2.0000 | 2.0104 | [1.9073, 2.1128] |
| 2 | Weight pi | 0.4000 | 0.3888 | [0.3673, 0.4102] |
| 2 | Shift kappa | 10.0000 | 10.0000 | [10.0000, 10.0000] |
| 2 | NB mean mu | 5.0000 | 5.2293 | [4.9679, 5.4989] |
| 2 | Shape r | 3.0000 | 3.3622 | [2.8346, 3.9984] |
| 2 | rtilde | 0.2500 | 0.2303 | [0.2001, 0.2608] |
| 2 | Success probability p | 0.3750 | 0.3907 | [0.3497, 0.4351] |
| 2 | Component mean kappa + mu | 15.0000 | 15.2293 | [14.9679, 15.4989] |
| Extra 1 | Weight pi | 0.0000 | 0.0004 | [7.170e-168, 0.0038] |
| Extra 1 | Shift kappa | N/A | 15.1136 | [0.0000, 30.0000] |
| Extra 1 | NB mean mu | N/A | 15.5826 | [0.7619, 30.2043] |
| Extra 1 | Shape r | N/A | 38.5772 | [1.0504, 76.0800] |
| Extra 1 | rtilde | N/A | 0.2486 | [0.0130, 0.4877] |
| Extra 1 | Success probability p | N/A | 0.3045 | [0.0435, 0.9296] |
| Extra 1 | Component mean kappa + mu | N/A | 30.6962 | [6.9221, 55.5731] |
| Extra 2 | Weight pi | 0.0000 | 1.977e-06 | [0.0000, 1.612e-06] |
| Extra 2 | Shift kappa | N/A | 15.3114 | [0.0000, 31.0000] |
| Extra 2 | NB mean mu | N/A | 15.6788 | [0.6882, 30.2669] |
| Extra 2 | Shape r | N/A | 17.6568 | [1.0567, 85.2917] |
| Extra 2 | rtilde | N/A | 0.2478 | [0.0116, 0.4862] |
| Extra 2 | Success probability p | N/A | 0.3048 | [0.0429, 0.9445] |
| Extra 2 | Component mean kappa + mu | N/A | 30.9902 | [6.0972, 55.4442] |
| Extra 3 | Weight pi | 0.0000 | 5.741e-08 | [0.0000, 1.838e-11] |
| Extra 3 | Shift kappa | N/A | 15.4236 | [0.0000, 31.0000] |
| Extra 3 | NB mean mu | N/A | 15.6713 | [0.8428, 30.2639] |
| Extra 3 | Shape r | N/A | 17.1730 | [1.0504, 68.1987] |
| Extra 3 | rtilde | N/A | 0.2512 | [0.0145, 0.4877] |
| Extra 3 | Success probability p | N/A | 0.2967 | [0.0432, 0.9189] |
| Extra 3 | Component mean kappa + mu | N/A | 31.0949 | [6.6442, 55.5204] |
| Extra 4 | Weight pi | 0.0000 | 1.123e-09 | [0.0000, 7.864e-21] |
| Extra 4 | Shift kappa | N/A | 15.9920 | [0.0000, 31.0000] |
| Extra 4 | NB mean mu | N/A | 15.4908 | [0.7975, 30.2513] |
| Extra 4 | Shape r | N/A | 17.1726 | [1.0533, 83.3500] |
| Extra 4 | rtilde | N/A | 0.2509 | [0.0119, 0.4870] |
| Extra 4 | Success probability p | N/A | 0.3012 | [0.0430, 0.9385] |
| Extra 4 | Component mean kappa + mu | N/A | 31.4828 | [6.7002, 55.5906] |

Combined weight of the four extras: posterior mean 0.0004; 95% credible interval [7.170e-168, 0.0038]. Their individual parameter summaries describe weight-ranked slots, not persistent component identities.

## Reproduction and checks

Run from the repository root:

```r
source("JOE revision/run-designs.R")
run_designs()
```

To regenerate only the report and plots from saved chains, without running any new chains:

```r
write_report(readRDS("JOE revision/results.rds"), "JOE revision")
```

The runner saves observations, full chains, retained draws, summaries, and R session details in [results.rds](results.rds). Seeds and settings are recorded above. It checks finite draws, normalized weights, PMF non-negativity and mass, and agreement between stored and independently recomputed log-likelihoods at ten retained draws per case. These numerical checks do not certify convergence.

Report generated: 2026-11-06. Sampling session: R version 4.5.3 (2026-03-11 ucrt). BayesMultiMode 0.7.1 loaded from the local working tree.
