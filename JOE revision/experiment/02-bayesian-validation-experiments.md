# Bayesian estimation: a one-hour practical check

Check whether the existing sampler mixes, gives similar answers across chains, and recovers a known distribution, its parameters, and its modes. This is a practical diagnostic study, not a calibration campaign or a comparison of new algorithms. No new sampler is required.

V1-V4 are experiment labels. Write components as $(\kappa,\mu,r)$: shift, unshifted NB mean, and shape. The total component mean is $\kappa+\mu$. $K$ is the fitted component limit; $L$ is the number of modes. R1-R3 identify reviewers.

## Experiments

| Experiment | Simulated truth and fitted model | What to check | Reviewer comments |
|---|---|---|---|
| V1. Single-component recovery | $n=1000$ from $(5,8,3)$; fit $K=1$ using the existing underlying sampler. | Recover shift, mean, shape, and the true distribution and modes. Compare posterior estimates and 95% intervals with truth. This is a basic recovery check, not an independent calculation of the posterior. | [R1 4: simulation evidence](../review/Reply-to-referee-1.tex#L222); [R2 5: validation](../review/Reply-to-referee-2.tex#L244). |
| V2. Clearly separated peaks | $n=1000$; weights $(0.6,0.4)$; components $(0,2,5)$ and $(12,3,3)$; fit $K=2$. | Check mixture-level mixing and agreement across chains. Recover weights and component parameters, the true mode count, and locations. Report posterior probability of the true count and location probabilities conditional on that count. | [R1 5: mode locations](../review/Reply-to-referee-1.tex#L242); [R2 5: mode recovery](../review/Reply-to-referee-2.tex#L244); [R2 7: convergence](../review/Reply-to-referee-2.tex#L257). |
| V3. Components do not equal modes | $n=1000$; equal weights; components $(0,5,4)$ and $(1,5,4)$; fit $K=2$. | Calculate true modes directly from the population probabilities. Check whether the fitted mixture recovers one peak without confident false multimodality. Assess mixture recovery; do not demand precise recovery of weakly identified component parameters. | [R2 1: mixtures versus multimodality](../review/Reply-to-referee-2.tex#L168); [R1 4: simulation properties](../review/Reply-to-referee-1.tex#L222). |
| V4. Does overfitting spoil inference? | Reuse exactly the V2 observations; fit $K=6$. | Compare mixture fit, mode probabilities, and chain agreement with $K=2$. Report extra-component weight without interpreting individual empty components. This is a small capacity check, not the referee's full five-component $K=5,10,15$ study. | [R3 2: component-limit sensitivity](../review/Reply-to-referee-3.tex#L161); [R3 3: mixture inference](../review/Reply-to-referee-3.tex#L165). |

## Settings and time budget

Use two independently seeded datasets per truth and four independent chain seeds per fit: 32 chains total. Start with 5,000 iterations per chain, discarding 2,500 warm-up draws. Use existing default priors and warm-up adaptation. Record settings, seeds, actual starting states, and runtime. The current sampler's starts follow its initialization policy; different seeds do not guarantee dispersed starts.

Timing check on this machine: 500 iterations took 0.50, 0.86, 0.56, and 1.43 seconds for V1-V4. Linear extrapolation gives about 5 minutes for the planned sampling, excluding analysis and reruns; it is not a convergence result. Reserve 40 minutes for sampling, 15 for summaries/plots, and 5 as a buffer. Use a 55-minute stop/checkpoint limit and report unfinished fits. Do not shorten chains to manufacture convergence.

## Mixing and convergence

Use trace plots, autocorrelation, rank-normalized $\widehat R$, and bulk/tail effective sample size (ESS) for log likelihood, mixture mean/variance, zero probability, and selected population-probability values. Aim for $\widehat R<1.01$ and ESS above 400 across chains; report failures. Check mode-count indicators and chain-specific mode probabilities too. Constant indicators have no informative R-hat; agreement alone does not establish mixing. Component labels are irrelevant to these diagnostics.

## Recovery and output

Show true versus posterior-mean probability curves with 95% pointwise bands; report total variation and omitted tail mass. Find modes for each retained draw on a checked integer range beyond all component modes. Include every component and count a flat peak once. Distinguish posterior mode probabilities from modes of the averaged distribution.

For V1/V2, show true values, posterior means, and 95% intervals for weights, shifts, means, and shapes. Transform $r=1/\tilde r-1$ in each draw before summarizing. Order distinguishable components by total mean only for reporting; this is not a label-switching investigation.

Deliver one diagnostic table, one recovery table, and trace/distribution plots. Two datasets per truth provide a useful check, not reliable repeated-sampling coverage or proof of convergence. Report poor recovery and unresolved mixing honestly.

Use the existing [sampler](../../R/gibbs_sfm_algos.R#L631); numerical checks of weights, support, and stored likelihoods remain required.