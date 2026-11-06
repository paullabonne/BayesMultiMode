# When shifted negative binomial mixtures help

Find out when shifted negative binomial (SNB) mixtures improve distribution estimation and mode inference, and when competitors do better. These proposed experiments test specific claims; they do not assume superiority. Existing designs are optional and should be used only if they answer the question in the table.

Write $\mathrm{SNB}(\kappa,\mu,r)=\kappa+\mathrm{NB}(\mu,r)$, where $\mu$ is the unshifted mean and $r$ controls dispersion. The total mean is $\kappa+\mu$; $r=\infty$ means Poisson. R1-R3 identify reviewers, with links to the relevant comments. All data are synthetic.

## Experiments

| Experiment | What to simulate | What would answer the question? | Reviewer comments |
|---|---|---|---|
| C1. Does the shift help? | Compare $\mathrm{SNB}(0,20,3)$ with $\mathrm{SNB}(15,5,3)$. Both have mean 20; the second concentrates above a nonzero starting point and has variance below its mean. | Compare shifted and unshifted NB mixtures. Does allowing shifts improve fit and mode locations? The unshifted truth checks whether extra flexibility has a cost. | [R1 4a: evidence for shifts](../review/Reply-to-referee-1.tex#L224); [R3 4: test the contributions](../review/Reply-to-referee-3.tex#L169). |
| C2. Does NB dispersion help? | Two equal-weight components $(10,5,100)$ and $(30,10,2)$. Repeat with both shapes 100, then with both infinite. | Compare SNB with shifted Poisson mixtures. Does SNB fit the broader tail without inventing extra modes? Does the advantage shrink for near-Poisson or Poisson truth? | [R1 4c: shifted Poisson comparison](../review/Reply-to-referee-1.tex#L230); [R3 4: evidence beyond fitted pictures](../review/Reply-to-referee-3.tex#L169). |
| C3. Can it detect weak or overlapping peaks? | Two components with shifts $(0,d)$, $d=3,10$; means $(2,5)$; shapes $(5,3)$; weights $(0.6,0.4)$ or $(0.95,0.05)$. Include a single $\mathrm{NB}(10,2)$ as a unimodal control. | Measure false extra modes, missed modes, and location uncertainty. Compute the true modes from the mixture probabilities; never assume each component produces a mode. | [R2 1: mixture versus multimodality](../review/Reply-to-referee-2.tex#L168); [R2 5: convincing synthetic evidence](../review/Reply-to-referee-2.tex#L244); [R1 5: mode locations](../review/Reply-to-referee-1.tex#L242). |
| C4. Can specialized zero models do better? | With probability 0.35 set the observation to zero. Otherwise draw from $\mathrm{NB}(8,2)$; in a second design, draw from that NB conditional on being positive. | Compare SNB with correctly specified zero-inflated and hurdle NB models. Assess zero probability, the positive tail, and modes. These are deliberate opportunities for competitors to win. | [R2 4: zero-inflated and hurdle alternatives](../review/Reply-to-referee-2.tex#L235); [R2 5: zero-inflated simulations](../review/Reply-to-referee-2.tex#L244). |
| C5. What if the component family is wrong? | Use $\mathrm{NB}(10,0.5)$, outside the implemented $r>1$ restriction. Separately use Poisson-lognormal counts: $\log\lambda\sim N(\log10-1/2,1)$ and $Y\mid\lambda\sim\mathrm{Poisson}(\lambda)$. | Compare with unrestricted-shape NB and Poisson-lognormal models. Report tail error, predictive accuracy, and false modes. This tests limits of the claimed flexibility, not sampler correctness. | [R2 1: scope of the model](../review/Reply-to-referee-2.tex#L168); [R3 1: alternative NB methods](../review/Reply-to-referee-3.tex#L157). |
| C6. Are modes caused by rounding? | Start with unimodal $X\sim\mathrm{Binomial}(40,0.25)$. Independently round each observation to the nearest multiple of five with probability $0,0.3,0.6$. | Compare with a binomial model that explicitly includes rounding. Evaluate the observed distribution separately from the underlying distribution. Recovering observed peaks alone cannot establish latent groups. | [R3 5: heaping and point masses](../review/Reply-to-referee-3.tex#L173). |

## Fair comparisons

Fit Poisson, shifted Poisson, unshifted NB, and SNB sparse mixtures throughout, plus each row's specialized competitors and a discrete kernel estimator. Use component limit $K=10$ initially. The [validation plan's V6 experiment](02-bayesian-validation-experiments.md#L16) covers $K=5,10,15$, prior sensitivity, and the Dirichlet process comparison; run it once for both studies.

Use $n=100,500,2000$, 500 independent datasets per setting, and 10,000 new test observations per dataset. Share training/test data across methods. Exclude 50 pilots used to set tuning and runtime budgets. Tune using training data only; match shared priors and compare other priors through the distributions they imply. Apply the validation plan's convergence checks and report failures and runtime.

Unshifted NB mixtures have variance at least as large as their mean. C1 tests a concrete limitation that shifts can overcome. NB tails are heavier than Poisson tails, but not power-law tails.

## What to report

Use two primary measures: predictive accuracy on new observations, measured by the mean log score, and accuracy of Bayesian mode-count probabilities, measured by Brier loss. The latter penalizes confident wrong answers. Average posterior-draw probabilities before taking predictive logs. For non-Bayesian estimators, report distribution and point-estimate errors separately rather than invent posterior probabilities.

Also report total variation, tail error, false/missed modes, location error, and interval coverage/width. Count flat peaks once; search beyond all component modes and check omitted tail mass. Keep all components. Do not treat the sample maximum as a population boundary. Any threshold for ignoring tiny modes needs separate, declared results.

Show paired differences with Monte Carlo uncertainty, adjusting prespecified multiple comparisons using Holm's method. Claim an advantage only where both primary outcomes support it. Summarize each experiment as: where SNB wins, where it ties or loses, and why. New comparison code is needed; no existing experiment is compulsory.