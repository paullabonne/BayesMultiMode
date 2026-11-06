# When shifted negative binomial mixtures help: a one-hour comparison

Compare the existing shifted negative binomial (SNB), unshifted NB, and shifted Poisson mixture estimators. Check distribution fit, predictive accuracy, and recovery of mode counts and locations. No new sampler, specialized competitor, or large replication grid is required.

Write $\mathrm{SNB}(\kappa,\mu,r)=\kappa+\mathrm{NB}(\mu,r)$: shift, unshifted mean, and shape. Total component mean is $\kappa+\mu$. C1-C6 are experiment labels; R1-R3 identify reviewers. These are proposed synthetic examples, not completed evidence of superiority.

## Experiments

| Experiment | One population setting | What to compare | Reviewer comments |
|---|---|---|---|
| C1. Nonzero shift | $\mathrm{SNB}(15,5,3)$. Mean 20, variance approximately 13.3. | Does SNB improve distribution fit and mode locations versus unshifted NB? Unshifted NB mixtures cannot have variance below their total mean. | [R1 4a: evidence for shifts](../review/Reply-to-referee-1.tex#L224); [R3 4: test the contributions](../review/Reply-to-referee-3.tex#L169). |
| C2. Unshifted control | $\mathrm{SNB}(0,20,3)$. Same mean as C1, but no shift. | Does the extra shift flexibility help, tie, or hurt when it is unnecessary? Do not assume SNB should win. | [R1 4a: shift comparison](../review/Reply-to-referee-1.tex#L224); [R1 4c: different generating processes](../review/Reply-to-referee-1.tex#L230). |
| C3. Narrow peak and broad tail | Equal weights; components $(10,5,100)$ and $(30,10,2)$. | Does SNB capture the broad tail better than shifted Poisson without inventing extra modes? Report both predictive fit and mode recovery. | [R1 4c: shifted Poisson comparison](../review/Reply-to-referee-1.tex#L230); [R3 4: numerical evidence](../review/Reply-to-referee-3.tex#L169). |
| C4. Poisson control | Equal weights; $10+\mathrm{Poisson}(5)$ and $30+\mathrm{Poisson}(10)$. | Does the dispersion advantage disappear under Poisson truth? Check whether SNB ties or loses to the simpler model. | [R1 4c: different generating processes](../review/Reply-to-referee-1.tex#L230). |
| C5. Overlapping components | Equal weights; components $(0,5,4)$ and $(1,5,4)$. | The population has one peak. Compare false extra modes and uncertainty about locations; two components must not automatically imply two modes. | [R2 1: mixtures versus multimodality](../review/Reply-to-referee-2.tex#L168); [R1 5: mode locations](../review/Reply-to-referee-1.tex#L242). |
| C6. Excess zeros | With probability 0.35 set $Y=0$; otherwise draw $\mathrm{NB}(8,2)$. | Compare zero probability, tail fit, and recovery of the zero and positive peaks using the three existing mixtures. This does not compare against a dedicated zero-inflated model. | [R2 5: zero-inflated simulations](../review/Reply-to-referee-2.tex#L244). |

## Settings and time budget

Use two independently seeded training datasets of $n=500$ per setting. Fit all three models to the same observations with $K=6$, four chain seeds, and 5,000 iterations per chain, discarding 2,500 warm-up draws: 144 chains total. Generate 5,000 independent test observations per dataset and share them across models. Use existing default priors; record differences between models, seeds, and actual starting states. These compare the implemented procedures, not identical posterior targets.

Short timing checks on all six settings took 0.42-2.51 seconds per 500-iteration fit across the three models. Extrapolation gives about 20 minutes for sampling, excluding analysis and reruns; these checks do not establish convergence. Reserve 40 minutes for sampling, 15 for analysis, and 5 as a buffer; stop/checkpoint at 55 minutes and report unfinished fits. Do not shorten chains merely to claim convergence.

Apply the [validation plan's mixture-level diagnostics](02-bayesian-validation-experiments.md#L1): traces, autocorrelation, rank-normalized $\widehat R$, and effective sample size for mixture quantities. Compare chain-specific mode probabilities. Label switching is not a diagnostic target. Report poor mixing and unreliable comparisons.

## Output

For each setting, show true versus estimated probability curves with 95% pointwise bands and a compact table: total variation, held-out mean log score, tail-probability error, posterior probability of the true mode count, and mode-location probabilities conditional on that count. Average posterior probabilities before taking predictive logs. Tabulate repeated test counts to avoid redundant probability evaluations.

Compute true modes from population probabilities and estimated modes for every retained draw. Include all components, count flat peaks once, and search beyond every component mode with checked tail mass. Do not use the observed maximum as a population boundary.

Report each dataset's paired model differences and Monte Carlo error, not significance tests or repeated-sampling coverage claims based on two datasets. End each example with where SNB wins, ties, or loses. Better fit alone does not guarantee better mode recovery.

Dedicated zero-inflated/hurdle models, heaping, alternative tail families, and broader prior sensitivity remain deferred referee requests. This comparison has its own one-hour budget, separate from validation.