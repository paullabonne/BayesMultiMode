# Bayesian estimation: a one-hour practical check

Check whether the existing sampler mixes, gives similar answers across chains, and recovers a known distribution, its parameters, and its modes. This is a practical diagnostic study, not a calibration campaign or a comparison of new algorithms. No new sampler is required for V1-V5 or V7.

V1-V7 are experiment labels. Write components as $(\kappa,\mu,r)$: shift, unshifted NB mean, and shape. The total component mean is $\kappa+\mu$. $K$ is the fitted component limit; $L$ is the number of modes. R1-R3 identify reviewers.

## Experiments

| Experiment | Simulated truth and fitted model | What to check | Reviewer comments |
|---|---|---|---|
| V1. Single-component recovery | $n=1000$ from $(5,8,3)$; fit $K=1$ using the existing underlying sampler. | Recover shift, mean, shape, and the true distribution and modes. Compare posterior estimates and 95% intervals with truth. This is a basic recovery check, not an independent calculation of the posterior. | [R1 4: simulation evidence](../review/Reply-to-referee-1.tex#L222); [R2 5: validation](../review/Reply-to-referee-2.tex#L244). |
| V2. Clearly separated peaks | $n=1000$; weights $(0.6,0.4)$; components $(0,2,5)$ and $(12,3,3)$; fit $K=2$. | Check mixture-level mixing and agreement across chains. Recover weights and component parameters, the true mode count, and locations. Report posterior probability of the true count and location probabilities conditional on that count. | [R1 5: mode locations](../review/Reply-to-referee-1.tex#L242); [R2 5: mode recovery](../review/Reply-to-referee-2.tex#L244); [R2 7: convergence](../review/Reply-to-referee-2.tex#L257). |
| V3. Components do not equal modes | $n=1000$; equal weights; components $(0,5,4)$ and $(1,5,4)$; fit $K=2$. | Calculate true modes directly from the population probabilities. Check whether the fitted mixture recovers one peak without confident false multimodality. Assess mixture recovery; do not demand precise recovery of weakly identified component parameters. | [R2 1: mixtures versus multimodality](../review/Reply-to-referee-2.tex#L168); [R1 4: simulation properties](../review/Reply-to-referee-1.tex#L222). |
| V4. Does overfitting spoil inference? | Reuse exactly the V2 observations; fit $K=6$. | Compare mixture fit, mode probabilities, and chain agreement with $K=2$. Report extra-component weight without interpreting individual empty components. This is a small capacity check; V5 adds the referee's five-component $K=5,10,15$ study. | [R3 2: component-limit sensitivity](../review/Reply-to-referee-3.tex#L161); [R3 3: mixture inference](../review/Reply-to-referee-3.tex#L165). |
| V5. Five-component truth | Five equally weighted SNB components with $\mu=3$, $r=5$. Use shifts $(0,12,24,36,48)$ for separated peaks and $(0,1,2,3,4)$ for overlap. Generate one dataset of $n=500$ per truth and fit each with $K=5,10,15$. | Compare mode-count probabilities and conditional locations across $K$, with mixing checks. Compute the true modes from the population PMF; five generating components do not imply five modes. This is a small robustness illustration, not a repeated-sampling study. | [R3 2: requested component limits](../review/Reply-to-referee-3.tex#L161); [R2 1: components versus modes](../review/Reply-to-referee-2.tex#L168). |
| V6. Prior sensitivity of mode inference | Reuse the two V5 datasets and their $K=10$ default-prior fits. Fit only the two alternative sparsity priors below to the same observations. | Determine whether prior changes alter $P(L\mid y)$, the most probable count, or locations conditional on that count. No density-shape, TV, tail-fit, or parameter-recovery comparison is required for this experiment. | [R2 6: prior specification](../review/Reply-to-referee-2.tex#L251); [R3 4: prior assumptions](../review/Reply-to-referee-3.tex#L169). |
| V7. Sample size with multimodality and large $K$ | Use V5's five equally weighted, separated components: shifts $(0,12,24,36,48)$, $\mu=3$, $r=5$. Fit nested samples with $n=100,500,1000$ at fixed $K=15$, reusing V5's $n=500$ fit. | Compare $P(L=5\mid y)$, the full mode-count distribution, and conditional mode locations and uncertainty as observations increase. Report missed or extra modes and mixing; do not require every realized result to improve monotonically. | [R1 4: simulation properties](../review/Reply-to-referee-1.tex#L222); [R2 5: mode recovery](../review/Reply-to-referee-2.tex#L244); [R3 2: oversized mixtures](../review/Reply-to-referee-3.tex#L161). |

## Settings and time budget

For V1-V4, use two independently seeded datasets per truth and four independent chain seeds per fit: 32 chains total. Start with 5,000 iterations per chain, discarding 2,500 warm-up draws. Use existing default priors and warm-up adaptation. Record settings, seeds, actual starting states, and runtime. The current sampler's starts follow its initialization policy; different seeds do not guarantee dispersed starts.

Timing check on this machine: 500 iterations took 0.50, 0.86, 0.56, and 1.43 seconds for V1-V4. Linear extrapolation gives about 5 minutes for the planned sampling, excluding analysis and reruns; it is not a convergence result. For this round, use one shared 55-minute limit for validation, optional comparisons, and reporting, not an hour per document. Prioritize V5, then V7, then V6; allow up to 40 minutes for sampling, 10 for summaries, and 5 for checkpointing. Time the first full V5 chain at each $K$ and the first V7 chain at $n=1000$, and use the remaining budget to decide which optional work fits. Do not start a fit or summary that is projected to exceed the deadline; checkpoint by minute 55 and report unfinished or unconverged fits. Do not shorten chains to manufacture convergence.

## Mixing and convergence

Use trace plots, autocorrelation, rank-normalized $\widehat R$, and bulk/tail effective sample size (ESS) for log likelihood, mixture mean/variance, zero probability, and selected population-probability values. Aim for $\widehat R<1.01$ and ESS above 400 across chains; report failures. Check mode-count indicators and chain-specific mode probabilities too. Constant indicators have no informative R-hat; agreement alone does not establish mixing. Component labels are irrelevant to these diagnostics.

## Recovery and output

For V1-V4, show true versus posterior-mean probability curves with 95% pointwise bands; report total variation and omitted tail mass. V5 and V7 need only the paired mode-inference tables described below. Find modes for each retained draw on a checked integer range beyond all component modes. Include every component and count a flat peak once. Distinguish posterior mode probabilities from modes of the averaged distribution.

For V1/V2, show true values, posterior means, and 95% intervals for weights, shifts, means, and shapes. Transform $r=1/\tilde r-1$ in each draw before summarizing. Order distinguishable components by total mean only for reporting; this is not a label-switching investigation.

Deliver one diagnostic table, one recovery table, and trace/distribution plots. Two datasets per truth in V1-V4 provide a useful check, not reliable repeated-sampling coverage or proof of convergence. Report poor recovery and unresolved mixing honestly.

Use the existing [sampler](../../R/gibbs_sfm_algos.R#L631); numerical checks of weights, support, and stored likelihoods remain required.

## V5-V7 quick checks

Use four chains per fit, with 5,000 iterations and 2,500 warm-up draws. V5 requires $2\times3\times4=24$ chains; V6 adds $2\times2\times4=16$, reusing the eight $K=10$ baseline chains. V7 adds $2\times4=8$ chains for $n=100$ and $n=1000$, reusing V5's four $n=500$, $K=15$ chains: 48 new chains in total. Do not cross the sample-size check with other $K$ values or priors, or add a replication grid. Reuse saved V1-V4 checks when their settings and package revision match; reruns are optional within the same deadline. The package remains the reference implementation. These additions require later extensions to [run-validation.R](run-validation.R), not changes in this documentation task.

For V5, keep componentwise priors and the sparsity hyperprior fixed across $K=5,10,15$; total Dirichlet concentration therefore changes with $K$. Report one table per dataset with $K$, the preferred mode count, $P(L=L_{\mathrm{true}}\mid y)$, and conditional mode locations, with Monte Carlo errors. Show whether the conclusions change across all three values, including $K=10$ versus $K=15$. Report sensitivity if found; two examples cannot establish universal robustness.

For V6, vary one prior choice at a time at fixed $K=10$. Gamma parameters below are shape and rate.

| Prior choice | Baseline and alternatives |
|---|---|
| Sparsity strength | $\alpha\sim\mathrm{Gamma}(1,200)$; change the rate to 50 and 800. |

Report only changes in $P(L\mid y)$, the preferred count, and conditional locations for counts in the union of the baseline and variant 95% credible sets, with Monte Carlo errors. Retain convergence diagnostics, but do not compare density shape, TV, tail fit, or parameter recovery across priors. Further prior variants, zero inflation, and heaping are deferred; the saved Poisson-prior check need not be rerun in this round.

For V7, use the first 100 observations of V5's separated-peak dataset, retain that entire dataset for $n=500$, and append 500 independently generated observations from the same truth for $n=1000$. Fix seeds before fitting. Keep $K=15$, the package prior rules, and chain settings unchanged; record any data-dependent prior bounds. Report one row per sample size with the preferred mode count, $P(L=5\mid y)$, a 95% credible set for $L$, and conditional-location summaries with Monte Carlo errors. These nested samples illustrate the effect of sample size under the package procedure, not repeated-sampling performance.