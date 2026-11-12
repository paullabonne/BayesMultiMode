# Econometrics response plan

The three-document structure works: this plan maps the objections, the [validation experiments](02-bayesian-validation-experiments.md#L1) check the Bayesian procedure, and the [comparison experiments](03-shifted-nb-comparison-experiments.md#L1) identify when shifted negative binomial (SNB) mixtures help. Simulations provide evidence within specified settings, not universal proof. Select experiments for the reviewers' questions; existing designs are optional.

## Revision priorities

| Concern | Referee comments | Required response |
|---|---|---|
| What is new? | [Editor 1](../review/Reply-to-editor.tex#L157), [2](../review/Reply-to-editor.tex#L164); [R1 1](../review/Reply-to-referee-1.tex#L155), [6](../review/Reply-to-referee-1.tex#L249), [8](../review/Reply-to-referee-1.tex#L260); [R2 1](../review/Reply-to-referee-2.tex#L168), [2](../review/Reply-to-referee-2.tex#L185); [R3 1](../review/Reply-to-referee-3.tex#L157) | Credit existing shifts, NB mixtures, and sampling methods. Demonstrate better computation and useful uncertainty about discrete modes. Explain the relationship to regression methods without claiming this paper provides regression. |
| What are we estimating? | [R1 5](../review/Reply-to-referee-1.tex#L242); [R2 1](../review/Reply-to-referee-2.tex#L168), [3](../review/Reply-to-referee-2.tex#L205); [R3 3](../review/Reply-to-referee-3.tex#L165) | Distinguish fitted components $K$, occupied components $K_+$, and modes $L$. Explain mode-count probabilities and mode locations conditional on the count. Component labels do not change the fitted distribution. |
| Does the code match the mathematics? | [R1 3a](../review/Reply-to-referee-1.tex#L206), [3d](../review/Reply-to-referee-1.tex#L212), [minor comments](../review/Reply-to-referee-1.tex#L276) | Reconcile success-probability notation, Gamma rate/scale, parameter bounds, allocations, and empty-component updates. Explain the bounded conditional when the usual inverse-Gamma distribution is unavailable. V1 checks recovery, not the mathematical correctness of every update. |
| Does estimation work in practice? | [R1 3](../review/Reply-to-referee-1.tex#L204), [4b](../review/Reply-to-referee-1.tex#L227), [4d](../review/Reply-to-referee-1.tex#L233); [R2 3](../review/Reply-to-referee-2.tex#L205), [7](../review/Reply-to-referee-2.tex#L257); [R3 3](../review/Reply-to-referee-3.tex#L165), [4](../review/Reply-to-referee-3.tex#L169) | Describe the implemented updates. V1-V4 check mixture mixing, chain agreement, parameter recovery, and mode recovery; reuse matching saved runs. V7 checks five-mode truth at $n=100,500,1000$ with $K=15$. C12 retains the pilot-inspired parameterization and sampler-efficiency design but is deferred beyond the quick check. |
| Do priors drive the result? | [R2 6](../review/Reply-to-referee-2.tex#L251); [R3 2](../review/Reply-to-referee-3.tex#L161), [4](../review/Reply-to-referee-3.tex#L169) | State the actual priors. V4 provides a small component-limit check. Prioritize V5: two datasets of five-component truth with $K=5,10,15$. V6 reuses their $K=10$ fits for two alternative sparsity priors, reporting mode inference only. V5-V7 share one 55-minute budget with optional comparisons and reporting. A Dirichlet process benchmark remains deferred. |
| When does SNB outperform? | [R1 4a](../review/Reply-to-referee-1.tex#L224), [4c](../review/Reply-to-referee-1.tex#L230); [R2 4](../review/Reply-to-referee-2.tex#L235), [5](../review/Reply-to-referee-2.tex#L244); [R3 4](../review/Reply-to-referee-3.tex#L169) | C1-C5 compare SNB, unshifted NB, and shifted Poisson; reuse matching saved runs. Test shifts, dispersion, and overlapping peaks, including controls where SNB need not win. C8, C9, and C11 are optional single-dataset additions within the shared deadline; C7 and C10 are reserve designs. C6, C13, zero-inflated/hurdle models, and heaping are deferred. |
| What can the results establish? | [R1 9](../review/Reply-to-referee-1.tex#L265); [R3 5](../review/Reply-to-referee-3.tex#L173) | Prediction of new independent counts is not temporal forecasting. Observed peaks alone do not identify latent groups or economic equilibria. A heaping-aware model remains a separate referee request, outside these compact studies. |

## First: clarify the model

Use $E(Y\mid j)=\kappa_j+\mu_j$ and $\operatorname{Var}(Y\mid j)=\mu_j+\mu_j^2/r_j$. Explain that $\tilde r\sim U(0,1/2)$ implies $p(r)=2/(1+r)^2$, $r>1$: proper does not mean uninformative.

Provide a short approximation argument: as $\mu\to0$, an SNB component approaches a point mass at its shift; enough components can approximate finite-support count laws. State the required capacity and support. This does not establish posterior consistency or reliable mode recovery near plateaus.

Use all components when computing each draw's distribution, including empty components with positive weights. Removing them changes what is estimated. Define flat peaks and boundary modes; finding modes after averaging distributions is not the same as calculating mode probabilities across draws. Probabilities that individual locations are modes need not sum to one.

## Then: validate and compare

Show that the implemented updates preserve the intended posterior before benchmarking them. Explain the joint shift-mean transformation and its prior and bounds. Sparse finite mixtures keep a fixed component limit while occupancy varies; weights are not exactly zero. Compare this with reversible-jump sampling without claiming either method is universally better. Parallelize independent chains and post-processing, not sequential Gibbs iterations.

Run computational validation before the comparative experiments. Retain unfavorable results and narrow claims accordingly. Economic application revisions remain outside these documents; this econometrics plan alone does not resolve the editor's relevance concern.

The immediate priority is the one-hour estimation check. Formal calibration and alternative samplers are not prerequisites for that check; broader referee requests remain separate revision tasks.

## Source basis

[Manuscript](../review/BCHLVD_BayesDens.tex#L93), the linked referee comments, [sampler](../../R/gibbs_sfm_algos.R#L631), and [mode inference](../../R/bayes_mode.R#L154).