# When shifted negative binomial mixtures help: a one-hour comparison

Compare the existing shifted negative binomial (SNB), unshifted NB, and shifted Poisson mixture estimators. Check distribution fit, predictive accuracy, and recovery of mode counts and locations. Reuse the existing C1-C5 results; additional comparisons are optional within the validation plan's shared 55-minute budget. C12-C13 are deferred.

Write $\mathrm{SNB}(\kappa,\mu,r)=\kappa+\mathrm{NB}(\mu,r)$: shift, unshifted mean, and shape. Total component mean is $\kappa+\mu$. C1-C13 are experiment labels; R1-R3 identify reviewers. These are proposed synthetic examples, not completed evidence of superiority.

## Experiments

| Experiment | One population setting | What to compare | Reviewer comments |
|---|---|---|---|
| C1. Nonzero shift | $\mathrm{SNB}(15,5,3)$. Mean 20, variance approximately 13.3. | Does SNB improve distribution fit and mode locations versus unshifted NB? Unshifted NB mixtures cannot have variance below their total mean. | [R1 4a: evidence for shifts](../review/Reply-to-referee-1.tex#L224); [R3 4: test the contributions](../review/Reply-to-referee-3.tex#L169). |
| C2. Unshifted control | $\mathrm{SNB}(0,20,3)$. Same mean as C1, but no shift. | Does the extra shift flexibility help, tie, or hurt when it is unnecessary? Do not assume SNB should win. | [R1 4a: shift comparison](../review/Reply-to-referee-1.tex#L224); [R1 4c: different generating processes](../review/Reply-to-referee-1.tex#L230). |
| C3. Narrow peak and broad tail | Equal weights; components $(10,5,100)$ and $(30,10,2)$. | Does SNB capture the broad tail better than shifted Poisson without inventing extra modes? Report both predictive fit and mode recovery. | [R1 4c: shifted Poisson comparison](../review/Reply-to-referee-1.tex#L230); [R3 4: numerical evidence](../review/Reply-to-referee-3.tex#L169). |
| C4. Poisson control | Equal weights; $10+\mathrm{Poisson}(5)$ and $30+\mathrm{Poisson}(10)$. | Does the dispersion advantage disappear under Poisson truth? Check whether SNB ties or loses to the simpler model. | [R1 4c: different generating processes](../review/Reply-to-referee-1.tex#L230). |
| C5. Overlapping components | Equal weights; components $(0,5,4)$ and $(1,5,4)$. | The population has one peak. Compare false extra modes and uncertainty about locations; two components must not automatically imply two modes. | [R2 1: mixtures versus multimodality](../review/Reply-to-referee-2.tex#L168); [R1 5: mode locations](../review/Reply-to-referee-1.tex#L242). |
| C6. Excess zeros (deferred) | With probability 0.35 set $Y=0$; otherwise draw $\mathrm{NB}(8,2)$. | Compare zero probability, tail fit, and recovery of the zero and positive peaks using the three existing mixtures. This does not compare against a dedicated zero-inflated model. | [R2 5: zero-inflated simulations](../review/Reply-to-referee-2.tex#L244). |
| C7. Separated unshifted components (pilot G3) | Weights $(0.6,0.4)$; components $(0,3,2)$ and $(0,30,5)$. | Can shifted Poisson recover broad peaks without extra modes when shifts are unnecessary? Compare against both NB models. | [R1 4c: multiple generating processes](../review/Reply-to-referee-1.tex#L230). |
| C8. Narrow shifted peak (pilot G5) | Equal weights; components $(0,15,10)$ and $(50,1,2)$. | Separate the need for a shift from the benefit of NB dispersion. Check recovery of the narrow peak near 50 as well as overall fit. | [R1 4a: evidence for shifts](../review/Reply-to-referee-1.tex#L224); [R1 4c: shifted Poisson comparison](../review/Reply-to-referee-1.tex#L230). |
| C9. Three nearby peaks (pilot G6) | Weights $(0.4,0.3,0.3)$; components $(0,2,5)$, $(7,3,5)$, and $(14,3,5)$. | Test recovery of several peaks and their conditional locations, including merging or inventing peaks. | [R1 4a](../review/Reply-to-referee-1.tex#L224); [R2 5: mode recovery](../review/Reply-to-referee-2.tex#L244). |
| C10. Three shifted-Poisson components (pilot P1) | Weights $(0.4,0.3,0.3)$; shifts $(0,7,14)$; each unshifted rate is 2.5. | Match the three-component structure of C9 under Poisson truth. Quantify the cost of unnecessary dispersion flexibility; SNB need not win. | [R1 4c: Poisson controls](../review/Reply-to-referee-1.tex#L230). |
| C11. Discretized normal mixture (pilot G8) | Equal weights; draw $X$ from $N(10,1.5^2)$ or $N(25,2^2)$, then $Y=\max(0,\lfloor X+0.5\rfloor)$. | Test both families outside their generating model. Keep settings where shifted Poisson wins. This is discretization to integer counts, not reporting heaping at multiples of five. | [R1 4c: different generating processes](../review/Reply-to-referee-1.tex#L230); [R3 4: broader evidence](../review/Reply-to-referee-3.tex#L169). |

C7-C11 draw on the [pilot designs](../Pilot_R1/code/s2_functions.R); the package remains the reference implementation. C6, zero inflation, and heaping are deferred. For C11, compute the true PMF with normal-CDF differences, pooling $X<0.5$ at zero to match the simulator.

## Settings and time budget

The C1-C5 protocol uses two independently seeded training datasets of $n=500$ per setting. Fit all three models to the same observations with $K=6$, four chain seeds, and 5,000 iterations per chain, discarding 2,500 warm-up draws: 120 chains total. Generate 5,000 independent test observations per dataset and share them across models. Use existing default priors; record differences between models, seeds, and actual starting states. These compare the implemented procedures, not identical posterior targets. Reuse saved results when their settings and package revision match; a full rerun is not part of this quick check.

Short timing checks on all six settings took 0.42-2.51 seconds per 500-iteration fit across the three models. Extrapolation gives about 20 minutes for sampling, excluding analysis and reruns; these checks do not establish convergence. This round shares the [validation plan's 55-minute limit](02-bayesian-validation-experiments.md): run V5 first, V7 next, then V6, and optional comparisons only if timed estimates leave room for reporting and checkpointing. Do not shorten chains merely to claim convergence.

Apply the [validation plan's mixture-level diagnostics](02-bayesian-validation-experiments.md#L1): traces, autocorrelation, rank-normalized $\widehat R$, and effective sample size for mixture quantities. Compare chain-specific mode probabilities. Label switching is not a diagnostic target. Report poor mixing and unreliable comparisons.

The commands below describe the existing C1-C6 runner, which still includes C6. It needs later amendments to omit C6 and run the additional experiments; no code changes or runs are part of this plan update.

For a selective rerun, archive the existing outputs, set `BAYESMULTIMODE_COMPARISON_MODEL`
to `SPois`, `SNB`, or `NB`, and run [run-comparison.R](run-comparison.R). The runner
reuses saved fits for the other models, checks that their test observations match,
and refreshes the tables and plots. Unset the variable to fit all three models.

For the shifted-Poisson prior-sensitivity rerun, also set
`BAYESMULTIMODE_COMPARISON_SPOIS_PRIOR=broader`. This uses a Gamma prior with
shape 1 and rate 0.05 for the Poisson rate, instead of shape 5 and rate 4.
Keep the same data, seeds, component limit, and chain lengths. Only the 48
shifted-Poisson chains are rerun; the SNB and NB results are reused. The figure
titles and chain-table columns `prior_l0` and `prior_L0` record the new prior.
This changes the comparison experiment, not the package defaults. Unset the
prior variable, or set it to `default`, to use the original prior.

For the C2 sample-size check, set `BAYESMULTIMODE_COMPARISON_SETTING=C2`,
`BAYESMULTIMODE_COMPARISON_N_TRAIN=5000`,
`BAYESMULTIMODE_COMPARISON_MODEL=SPois`, and
`BAYESMULTIMODE_COMPARISON_SPOIS_PRIOR=broader`. Set
`BAYESMULTIMODE_COMPARISON_OUT_DIR` to
`JOE revision/experiment/comparison-results-c2-n5000-spois` before running
[run-comparison.R](run-comparison.R). This fits eight chains across two
datasets without reusing the 500-observation fits. The test samples, seeds,
component limit, iterations, and burn-in remain unchanged. A separate output
directory is required when changing the training sample size; it contains
only the selected model's results and figures. The original comparison
outputs remain intact.

## Output

For each setting, show true versus estimated probability curves with 95% pointwise bands and a compact table: total variation, held-out mean log score, tail-probability error, posterior probability of the true mode count, and mode-location probabilities conditional on that count. Average posterior probabilities before taking predictive logs. Tabulate repeated test counts to avoid redundant probability evaluations.

Compute true modes from population probabilities and estimated modes for every retained draw. Include all components, count flat peaks once, and search beyond every component mode with checked tail mass. Do not use the observed maximum as a population boundary.

Report each dataset's paired model differences and Monte Carlo error, not significance tests or repeated-sampling coverage claims based on two datasets. End each example with where SNB wins, ties, or loses. Better fit alone does not guarantee better mode recovery.

Dedicated zero-inflated/hurdle models, heaping, and alternative tail families remain deferred referee requests. There is no separate additional hour for comparisons.

## Additional comparison settings

If time remains after V5-V7, run C8, then C9, then C11, with one dataset of $n=500$ per setting, $K=6$, and four chains per model: 12 chains per setting, at most 36. Use the C1-C5 chain lengths and paired test observations. Leave C7 and C10 as reserve designs; do not expand to more datasets or sample sizes. The [validation plan](02-bayesian-validation-experiments.md) covers $K=5,10,15$, the small mode-only prior check, and V7's $n=100,500,1000$ comparison at $K=15$. Report any comparisons left unfinished at the shared deadline.

Report each dataset's paired differences, mode-count probabilities, and conditional locations with Monte Carlo errors; do not estimate repeated-sampling success rates from these few examples. Show conditional locations for competing posterior counts, not just the true count. Use the package's full mixture throughout, including empty components with positive weights.

For C1, C8, and C9, add the pilot's analytic approximation check: unshifted NB mixtures satisfy $q(y)\leq\mathrm{Poisson}(y;y)$ for $y\geq1$, so TV error against truth $p$ is at least $\sum_{y\geq1}\max\{p(y)-\mathrm{Poisson}(y;y),0\}$. Report this lower bound alongside simulation errors; a zero bound is inconclusive.

## C12. Parameterization and sampler efficiency (deferred)

Use the [pilot benchmark design](../Pilot_R1/code/run_nbcomp.R) for [R1 4b/4d](../review/Reply-to-referee-1.tex#L227) and [R3 3/4](../review/Reply-to-referee-3.tex#L165). Start with one C9 dataset at $n=500$, $K=6$, and four chains per method. This requires alternative updates and is outside the current under-one-hour check.

| Method | Comparison |
|---|---|
| Package SNB sampler | Reference implementation. |
| Random walks in $(\mu,\tilde r)$ | New coordinates without the augmented mean update. |
| Random walks in $(p,r)$, $p=r/(r+\mu)$ | Traditional coordinates under the same transformed prior and bounds. |
| Zhou-style augmentation with MH correction | Alternative update targeting the package posterior. |
| Zhou-style augmentation under its own prior | Optional mode-inference sensitivity result, excluded from same-target efficiency rankings. |

Hold data, starting states, other update blocks, and post-processing common where mathematically valid; verify posterior targets and Jacobians before timing. Report acceptance, R-hat, ESS per second for label-invariant quantities, and Monte Carlo errors for mode probabilities. Time all methods consistently on the same hardware, including warm-up. Distinguish gains from parameterization and augmentation; do not assume the package beats every alternative. Alternative updates require later implementation, not replacement of the package with pilot code.

## C13. Histogram peaks under unimodal truth (deferred)

For [R1 6](../review/Reply-to-referee-1.tex#L249), use the [pilot histogram experiment](../Pilot_R1/code/pilot_hist.R): 10,000 samples from $\mathrm{NB}(\mu=5,r=5)$ at each of $n=270,411,4691$. With unit-width integer bins and plateaus counted once, report the frequency of two or more peaks and its Monte Carlo error. Report the pilot's 1%-mass rule separately as a histogram heuristic.

Fit package SNB with $K=6$ and four chains to the first two prespecified samples at each size. Compare the histogram and posterior preferred mode counts on these same samples, and report $P(L>1\mid y)$. The histogram's 1% rule does not apply to package mode inference. This study is outside the current under-one-hour check.