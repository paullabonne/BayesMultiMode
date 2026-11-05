# JOE Revision: Simulation Evidence and Work Plan

Status reviewed: 5 November 2026. Code baseline: `e0d836f`.

Scope revised after an independent review: use three methods and six scenarios in the core study, with targeted computational and sensitivity experiments rather than a full factorial design.

Phase 1 is closed as a pilot assessment; Phase 2 preparation, smoke testing, and calibration are implemented and executed. See [simulation/README.md](simulation/README.md) for the closeout and current results. The matched prior comparison did not support blaming all failures on sparsity. Phase 2 calibration passes 14/18 mean/zero-mass checks and 12/18 modal checks; the full production comparison has not started.

## Recommendation

Build a repeated simulation study around three separate questions: does allowing component-specific shifts improve count-distribution approximation; does negative binomial flexibility improve mode inference relative to shifted Poisson mixtures; and does the new parameterisation improve MCMC efficiency for the same posterior? Add a focused sensitivity study of the sparse finite mixture component cap and priors.

The strongest case is accurate reconstruction of the overall mixture distribution and uncertainty-aware inference about its population modes, with evidence showing when the method helps and when simpler methods suffice. A visually convincing fitted distribution or a faster sampler alone does not establish that case. Simulations can support the methodological claims; they cannot by themselves resolve the editor's request for economic relevance or establish novelty in the literature.

## Inference Target and Location

The target is the overall mixture PMF and its features, especially the number and locations of modes and predictive probabilities. It is not identification of the individual generating components, their parameters, their labels, or observation allocations. Different component decompositions that reconstruct the same mixture PMF and modes count as successful fits.

Do not add component-parameter recovery, true-component-count accuracy, classification accuracy, or relabelling studies. Component parameters and allocations are computational devices. Occupied-component counts may be saved as secondary sparsity/cost diagnostics, not scored against the generating component count. Any location matching below refers to mixture modes, never mixture components.

All new simulation scripts, study-specific helpers, settings, saved results, figures, and computational records belong in the separate `JOE revision/simulation/` folder. Keep this plan in its current location. Leave the existing [simulation.R](simulation.R) as an exploratory reference; develop the reproducible study in the new folder, reusing existing package functions and public generators without moving package code or sampling helpers. Make package loading and output paths explicit so the study does not depend on an interactive working directory.

## Where We Are

| Area | Evidence currently in the repository | Status |
| --- | --- | --- |
| Proposed model and main internal benchmarks | [../R/bayes_fit.R](../R/bayes_fit.R) supports unshifted and shifted negative binomial and Poisson mixtures. | Available for a controlled comparison. |
| Posterior mode inference | [../R/bayes_mode.R](../R/bayes_mode.R) and [../R/mix_mode.R](../R/mix_mode.R) provide the mode-inference machinery. Commit `79f2613` fixes discrete mode detection for unrelated equal-height points. | Independent truth, boundary, plateau, and label-invariance checks pass; pilot modal convergence remains insufficient in many groups. |
| Negative binomial sampler | Recent commits add a Gaussian proposal (`32a8dd9`) and fix mean bounds and empty-component updates (`77f5768`). [../tests/testthat/test-gibbs_sfm_algos.R](../tests/testthat/test-gibbs_sfm_algos.R) includes support, acceptance-ratio, burn-in tuning, and empty-component regression tests. | Implementation work has progressed; comparative efficiency has not been demonstrated. |
| Log-likelihood aggregation | Commit `e0d836f` fixes aggregation in all six samplers; the regression test checks a sum of observation-level log mixture probabilities. | Fixed in code. Recompute any affected results; this is not itself a predictive comparison. |
| Revision simulation | [simulation.R](simulation.R) generates one shifted-NB example and fits the proposed model. Its closing comments identify shift benefits, shifted Poisson comparison, and MCMC performance as work to do. | Exploratory, not a completed simulation study. |
| Phase 1 pilot | [simulation/run-phase1.R](simulation/run-phase1.R) ran three scenarios, four methods, ten datasets per scenario, and four chains per fit. [simulation/diagnose-phase1.R](simulation/diagnose-phase1.R) adds two longer-chain checks and mixture-level traces. | Execution complete; production remains blocked by convergence, not dependencies. |
| Manuscript evidence | [BCHLVD_BayesDens.tex](BCHLVD_BayesDens.tex) contains likelihood-surface illustrations and empirical NB-versus-Poisson comparisons. | These do not replace repeated experiments with known population modes. |
| Referee responses | [Reply-to-referee-1.tex](Reply-to-referee-1.tex), [Reply-to-referee-2.tex](Reply-to-referee-2.tex), and [Reply-to-referee-3.tex](Reply-to-referee-3.tex) retain simulation and efficiency placeholders. [Reply-to-editor.tex](Reply-to-editor.tex) has empty replies. | The requested evidence and the corresponding responses remain unfinished. |
| External fitting examples | [../external_comp.Rmd](../external_comp.Rmd) contains external mixture-fitting examples. | Useful integration examples, not evidence of an equivalent-posterior NB sampler comparison. |

The current simulation script cannot be treated as a clean-session experiment. It references `rtilde` and `bayesmix` before creating them, uses fixed component sample counts rather than random mixture allocations, and indexes `loglik[burnin:5000]` despite setting burn-in to 9000. It has no repeated fits of competitors, saved performance summaries, or Monte Carlo uncertainty. Preserve it as a reference and implement the reproducible study separately in `JOE revision/simulation/`.

The status table above records the original review baseline. The Phase 1 closeout below records subsequent implementation and runs; it does not certify convergence of all runs or reproduce the manuscript's empirical results.

### Current Implementation Constraints

- [../R/gibbs_sfm_algos.R](../R/gibbs_sfm_algos.R) restricts `rtilde` to `(0, 0.5)`, hence NB size `r > 1`. Correctly specified NB scenarios and the equivalent-posterior sampler comparison must respect this restriction. Cases with `r <= 1` belong in an explicitly labelled misspecification study unless extending the model is separately agreed.
- NB mean draws use a fixed full-data bound, and the public NB prior settings currently expose `a0`, `A0`, and proposal scale `sig`, not arbitrary mean/dispersion prior families. Do not promise a broad prior comparison without deciding what additional implementation it requires.
- The unshifted Poisson default rate prior uses `1.1 / median(data)`. A zero-heavy dataset with median zero cannot use that default. Specify valid, documented benchmark priors before the structural-zero comparison rather than treating a default-setting failure as evidence against the model.
- [../R/bayes_mode.R](../R/bayes_mode.R) defaults to the observed data range. Population-mode assessment needs an explicit range including zero and a checked upper tail; otherwise unobserved boundary or tail modes can be excluded by construction.
- [../R/bayes_info.R](../R/bayes_info.R) evaluates the training data and averages log likelihoods over retained draws. It is not the proposed held-out posterior predictive score. [../R/bayes_mixture.R](../R/bayes_mixture.R) retains the supplied full `loglik` vector while trimming `mcmc`; align retained-draw indices explicitly whenever using it.
- Phase 1 repairs empty shifted-Poisson prior refresh, shift likelihood arithmetic, support-respecting shift proposals, and concentration log ratios. Dedicated regression and finite-grid balance checks supplement the log-likelihood test, which mocks the shift and concentration updates. Count-allocation underflow/support failures now stop explicitly rather than triggering uniform fallback allocations.
- The intended concentration prior is `a0 = 1, A0 = 200`. Count samplers now use the manuscript's allocation-conditioned concentration target, integrating out weights and immediately redrawing them. This avoids failures from rounded-zero Dirichlet weights without changing the model or flooring weights. The original pilot's `a0 = 10, A0 = 50` results remain historical, not the default-method benchmark.
- Count samplers initialise with `kmeans(..., centers = K)`. Check datasets with fewer distinct observations than `K`, and all-zero datasets, explicitly. Do not silently regenerate failed datasets or interpret an initialisation failure as poor statistical performance.

The normal `Rscript` command works and all required imports are installed. Sampler, independent PMF/mode, quadrature, prior-integration, and update-order checks pass. Both one-component coordinate samplers agree with independent quadrature and pass R-hat/ESS checks. The fixed-complexity traditional mixture adapter also runs successfully.

## What the Referees Require

| Request | Evidence to deliver |
| --- | --- |
| R1, Comment 4(a): benefit of shifts | Shifted versus unshifted fits to the same simulated data, with distributional and modal errors. |
| R1, Comment 4(b,d); R3, Comments 3-4: parameterisation and efficiency | Same model and posterior under traditional and proposed coordinates; posterior agreement, effective sample size, and efficiency per second. |
| R1, Comment 4(c): shifted NB versus shifted Poisson | Several generating processes, repeated datasets, held-out predictive scores, and mode-count and location performance. |
| R2, Comment 5: zero inflation and multimodality | Known unimodal and multimodal populations, including structural-zero cases and misspecification. |
| R3, Comment 2: component cap | A five-component generating process fitted with caps of 5, 10, and 15; assess stability of the overall PMF and its modes, not recovery of five components. |
| R2, Comments 6-7; R3, Comment 4 | Targeted prior sensitivity and convergence/efficiency reporting. |
| R2, Comment 8; R3, Comment 5 | Explicit response on external zero-inflated/hurdle benchmarks and rounding/heaping in the retained applications; not an automatic requirement to implement every alternative. |
| Editor; R1, Comments 1, 6, and 9 | A separate, concrete explanation of what modal uncertainty changes in the financial/economic applications. |

## Study Design

### 1. Separate the Model Benefits

Use three methods in production: `neg_binomial`, `shifted_neg_binomial`, and `shifted_poisson`. Retain unshifted `poisson` in the small pilot as a control, not throughout the production grid. Use the same component cap and a documented concentration-prior convention.

Compare shifted NB with unshifted NB to assess the benefit of shifts, and shifted NB with shifted Poisson to assess dispersion flexibility. These contrasts do not estimate a full shift-by-dispersion interaction. Avoid claiming that unshifted mixtures cannot approximate a distribution in principle: the practical question is approximation and inference with finite data and a specified modelling budget.

Match Poisson rate priors in the pilot where possible, and inspect induced priors on observable means and PMFs across methods. Uniform `rtilde` induces a size prior proportional to `1/(1+r)^2` on `r > 1`, not a flat size prior. Where priors cannot be matched, describe results as model-plus-prior performance and use one targeted sensitivity check; do not attribute all differences solely to shifts or dispersion.

### 2. Use a Small, Informative Scenario Set

Start with the following families. Fix exact parameters before the production runs, after checking the population PMF and its actual mode set.

| Scenario | Purpose |
| --- | --- |
| Unshifted Poisson, one clear mode | Simple control: measure false multimodality and the cost of unnecessary flexibility. |
| Two overlapping, overdispersed NB components with one population mode and a long right tail | Test spurious tail modes and establish explicitly that the number of components is not the number of modes. |
| Shifted Poisson mixture with two clear modes | Control in which the main competitor is correctly specified. |
| Shifted NB mixture with two clear modes and overdispersion | Test the combined benefit of shifts and dispersion flexibility. |
| Same basic mixture with a weak or closely spaced second mode | Test detection limits and whether posterior uncertainty reflects the difficulty. |
| Structural point mass at zero plus an NB count distribution | Address zero inflation with a generating process beyond an ordinary NB mixture. |

Keep the five-component generating process separate, for the targeted cap experiment only. Start with one application-relevant sample size, provisionally 1000; do not cross every scenario with both 250 and 1000. Fix the size and exact parameters before production based on the application and pilot cost, not favourable comparative results.

Use economically relevant count scales and zero proportions where suitable, informed by the loan application. If retaining the inflation-expectations application, give an explicit response on rounding/heaping: peaks introduced by response rounding should not automatically be interpreted as distinct belief groups. Start with a limited application check or discussion, adding a dedicated simulation only if needed to support the retained claims. Call these application-motivated designs, not validated economic mechanisms.

An external structural-zero comparator needs an explicit author decision, not a second full benchmark programme. If included, start with one established zero-inflated NB implementation on the structural-zero scenario. Do not automatically fit both ZINB and hurdle NB: they can represent the same marginal PMF with corresponding parameters. If choosing a hurdle model instead, specify its positive-count distribution and include truncation normalisation when scoring. These methods are not currently provided by `bayes_fit()`; verify the chosen implementation on a simple case and report only the uncertainty measures it actually supplies.

### 3. Define Truth and Outcomes Before Fitting

- Define the population mode set from the generating PMF, not the generating component count or the sample histogram. State how zero-boundary modes, adjacent plateaus, and ties are counted.
- The primary modal estimand is the mode set of each posterior draw's PMF. Modes of the averaged predictive PMF are a different estimand; show them only as a separately labelled descriptive summary. Location inclusion probabilities need not sum to one across locations.
- Check truth independently of the mode-finding code being evaluated. For these generating families, extend the integer grid beyond every component's rightmost mode, where all components decrease, and check neighbouring points explicitly. Negligible omitted tail mass alone cannot rule out a tiny remote mode. Check omitted mass separately for distributional-error calculations and apply equivalent tail checks to fitted draws.
- Count an adjacent plateau as one mode and predeclare a location convention, such as its leftmost integer. Apply it consistently to truth, posterior draws, matching, and coverage. If introducing a practical prominence threshold, report it separately from all mathematical modes rather than silently discarding small peaks.
- Generate random mixture allocations, then observations, reusing suitable existing public generators. Fixed component counts may be a supplementary design but should be labelled as conditional on those counts.
- Use independent training and test samples. Compare posterior predictive distributions on the same test observations, with stable log-sum-exp arithmetic. Preserve genuine zero predictive probabilities as infinite log-score penalties; do not replace them with arbitrary floors. Distinguish these from numerical underflow.
- Primary outcomes: total variation distance between true and fitted mixture PMFs with tail handling stated; mean held-out log predictive density per observation; mode-count Brier score; repeated-sample coverage and size of a stated credible set for the mode count; and mode-location error with unmatched/missed modes reported. None measures recovery of individual components.
- Keep mode-count accuracy, false/missed multimodality, and mean probability assigned to the true count as secondary summaries. For methods without a posterior, report the outcomes they actually supply rather than inventing comparable credible probabilities.
- Report location error conditional on the correct count using ordered locations or a fixed matching rule, alongside the unconditional frequency of count/matching failures. Any location coverage analysis must state whether coverage is joint or marginal and conditional on the mode count. Fixed-truth repeated-sample coverage is an empirical performance measure, not guaranteed Bayesian calibration.
- Score the average predictive PMF, not the average log likelihood across posterior draws; label this as out-of-sample density prediction, not a time-series forecasting exercise.
- Keep posterior occupied-component counts separate from posterior mode counts. Use the former only to explain computational sparsity; do not treat agreement with the generating component count as evidence of mixture reconstruction.

### 4. Isolate the Parameterisation Benefit

Establish baseline feasibility in Phase 1, before production. Run a separate computational experiment on a single NB component and a small fixed-complexity mixture. Compare the proposed coordinates with a documented traditional NB sampler; the current external examples are not a ready-made NB baseline. If the equivalent-posterior baseline cannot be delivered, resolve the scope and retained computational claims with the authors before committing to a large study.

The two samplers must target the same posterior distribution on the same physical parameters. Transform priors and include the relevant Jacobian; independent uniform priors in different coordinates generally define different models. Hold shifts, component structure, data, and augmentation fixed wherever possible. If the baseline also changes augmentation or update strategy, describe the result as a comparison of algorithms rather than a pure coordinate effect.

Give both methods a reasonable tuning budget, freeze adaptation before retaining draws, and check posterior agreement before comparing speed. Include strongly overdispersed and near-Poisson cases, subject to the model's supported parameter domain. Do not add a fixed-versus-adapted proposal factorial by default; use it only as a small follow-up if the retained claim specifically concerns adaptation.

Use several independently initialised chains. For mixture fits, report rank-normalised split R-hat, bulk/tail ESS, Monte Carlo standard errors, and ESS per second for label-invariant mixture quantities: PMF values, overall mean, tail probabilities, mode count, and mode indicators. Occupied-component counts and acceptance rates are secondary computational diagnostics. Do not require recovery or matching of component parameters; unprocessed component traces are not sufficient evidence of poor mixing because labels may switch.

Record acceptance and proposal scales for the NB dispersion update, shifts, and concentration parameter as applicable. Audit what the samplers currently return: printed acceptance rates are not a complete saved diagnostics interface. State whether timing covers fitting only or also mode extraction. Keep detailed timing and tuning records with the revision experiments; include only the computational reporting needed for the paper or appendix.

### 5. Test Sparse-Mixture Robustness

For the five-component scenario, compare `K = 5, 10, 15` on 50 identical datasets using the proposed method. Compare overall PMF reconstruction, mode inference, predictive accuracy, diagnostics, and cost. Occupied components describe sparsity only; recovering the five generating components is not a success criterion. This is a targeted experiment, not another full method/scenario/sample-size grid.

Hold one concentration-prior convention fixed for the cap experiment and state whether it concerns `e0` or total concentration `K * e0`. Separately test two alternative shrinkage settings for the proposed method on 25 datasets in each of two preselected difficult scenarios, reusing the main-study baseline fits. Investigate the other cap/concentration convention only if needed to explain sensitivity. Do not promise configurable mean/dispersion priors without additional agreed implementation.

Report any cap or prior sensitivity. One five-component generating process cannot establish a universal cap rule; frame any recommendation as applying to the tested settings. DP, RJ-MCMC, HMC, and variational methods are suggested alternatives to address explicitly, not simultaneous required implementations. Defer them unless an author decision or a specific unresolved claim makes one necessary.

## Initial Run Budget

Assume four independently initialised chains per Bayesian fit. Chain length and elapsed cost remain to be determined by the pilot.

| Experiment | Design | Model fits | Chain runs |
| --- | --- | --- | --- |
| Pilot | 3 scenarios x 4 methods x 10 datasets x 1 sample size | 120 | 480 |
| Core comparison | 6 scenarios x 3 methods x 100 datasets x 1 sample size | 1800 | 7200 |
| Component cap | 1 five-component scenario x 3 caps x 50 datasets, proposed method only | 150 | 600 |
| Parameterisation | 2 dispersion regimes x 2 complexities x 5 datasets x 2 samplers | 40 | 160 |
| Shrinkage sensitivity | 2 alternative settings x 2 scenarios x 25 datasets, proposed method only | 100 | 400 |
| Initial total | Excludes correctness checks, retries, and external comparators | 2210 | 8840 |

The original eight-scenario, two-size, four-method grid with 200 datasets would require 51200 chain runs before add-ons. The revised budget retains the required questions without that full factorial cost. Pilot datasets are separate from production; sensitivity baseline fits are reused from production.

With 100 production replications, a proportion has worst-case Monte Carlo standard error 0.05; this cannot resolve small differences reliably. Predeclare precision targets and expansion rules before production, and add replications regardless of which method leads. A second size on only the weak-mode and structural-zero scenarios would add 600 fits / 2400 chains at 100 replications. Broader sizes, external comparators, adaptation contrasts, and additional algorithms require an explicit scope decision and revised budget.

## Ordered Work Plan

### Phase 1: Establish a Trustworthy Pilot

- [x] Resolve the local R execution failure: the normal `Rscript` command works. Focused sampler and independent discrete-mode checks pass; session information is saved.
- [x] Complete package-scoped sampler/discrete-mode tests and independent study checks with the installed imports. The package's flat-mode warning is expected for the plateau test.
- [x] Audit the NB prior/support restrictions, mean and shift bounds, returned diagnostics, and burn-in indexing against the manuscript. Record notation and shift-update discrepancies in the study README for the manuscript revision; the audit does not establish unrun sensitivity claims.
- [x] Resolve empty shifted-Poisson prior refresh and check shift transition balance on a small finite grid. Compare shift/concentration ratios with stable log-scale references; invalid count allocations now fail explicitly.
- [x] Check both one-component NB samplers against independent quadrature as a correctness check only. Both pass mixture-function R-hat/ESS and quadrature-agreement gates; boundary/plateau and low-distinct-count checks pass.
- [x] Establish a feasible traditional-coordinate baseline with the same target posterior and transformed prior. The fixed-complexity mixture adapter returns finite draws; this establishes feasibility, not a speed advantage.
- [x] Implement the clean-session runner in `JOE revision/simulation/`, with explicit settings, package/output paths, random allocations, saved failures, and no references to old fitted objects. Use [simulation.R](simulation.R) only as an exploratory reference; existing package code and helpers remain in place.
- [x] Implement and independently verify PMFs and true mode sets for a simple Poisson control, a one-mode NB tail case, and a shifted bimodal NB case. The tests also verify invariance to relabelling or splitting equivalent components.
- [x] Fit the four existing mixture families to ten datasets per pilot scenario at one size: all 480 chains and 120 group evaluations completed, with saved seeds, draws, scores, failures, and diagnostics. Only 29 groups passed the convergence gates; the rest were retained and reported.
- [x] Test longer chains on two selected failures without changing data, priors, seeds, or gates. At 12000 iterations / 6000 burn-in, the selected shifted-NB group passes (maximum R-hat 1.004, minimum bulk ESS 995); the shifted-Poisson group still fails (R-hat 1.257, ESS 11).
- [x] Run the intended-prior check without confounding sampler changes: 96 chains compare both priors using the corrected allocation-conditioned concentration update, the same data/seeds, and the original chain length. Both priors pass 4/12 mean/zero-mass checks; sparsity alone is not a general explanation.
- [x] Close the pilot assessment with separate PMF/modal diagnostics, between-chain PMF distances, modal-probability MCSE, and retained failures. Carry unresolved precision issues into Phase 2 calibration rather than describing a combined flag as failure of the whole mixture model.

Phase 1 is closed as an exploratory assessment with unresolved diagnostics explicitly handed to Phase 2 calibration; the original criterion that every pilot summary be convergence-validated has not been met and is not waived for production. Independent correctness and baseline feasibility are established. Discrete mode summaries require finite R-hat and bulk ESS; undefined discrete tail ESS and constant indicators are reported separately. Continuous summaries still require both bulk and tail ESS. No component recovery is assessed.

### Phase 2: Freeze and Run the Main Comparison

- [x] Implement the six-scenario protocol in [simulation/phase2.R](simulation/phase2.R): one sample size, three methods, intended concentration prior, explicit truth/metrics/precision settings, and failure-preserving paired comparisons. Hash settings and code before execution.
- [x] Generate all 600 training/test datasets with seeds separate from Phase 1. Verify each population PMF/mode set independently and generate the six-panel population figure.
- [x] Run the 36-chain execution smoke test and 72-chain calibration with four chains per method at 12000 iterations / 6000 burn-in. No fitting/scoring errors; 14/18 mean/zero-mass gates, 12/18 modal gates, and 11/18 combined gates pass.
- [ ] Resolve the seven combined calibration failures under the recorded protocol before the full run. PMF agreement and modal uncertainty are separate requirements; do not drop failed datasets or weaken thresholds to make them pass.
- [ ] Run 100 production replications per scenario initially. Record whether an external structural-zero comparator is included or deferred; do not add both ZINB and hurdle models by default.
- [ ] Finalise retained chain lengths or a validated local sampling improvement using Phase 2 calibration. The candidate 12000/6000 schedule improves diagnostics but is not sufficient for every method/scenario.
- [ ] Save scenario/replicate/method/chain seeds, training/test data or regeneration settings, fit settings, code commit, session information, and resumable summaries under `JOE revision/simulation/`. Record failed or unconverged fits; do not silently drop them.
- [ ] Report paired method differences with Monte Carlo standard errors, failure rates, and secondary false/missed-mode summaries. Expand replications only under the predeclared precision rules.
- [ ] Use a small selection of representative plots, not selected successes, to explain the repeated results. Add a second size only to the weak-mode and structural-zero cases if an agreed precision or scientific question requires it.

Completion criterion: tables establish where the proposed method improves overall mixture reconstruction and modal inference, where it ties or loses, and how certain those comparisons are. Recovery of individual components is not assessed.

### Phase 3: Complete the Computational and Sensitivity Evidence

- [ ] Run the already validated traditional-versus-proposed coordinate experiment over two dispersion regimes, two complexities, and five datasets, with saved acceptance, R-hat, ESS, and timing summaries.
- [ ] Run the five-component `K = 5, 10, 15` experiment on 50 datasets and the two-setting, two-scenario shrinkage sensitivity on 25 datasets each.
- [ ] Explain any cap sensitivity before adding more cap/prior combinations; restrict recommendations to the tested settings.
- [ ] Give an explicit rounding/heaping response if retaining the inflation application. Record author decisions on the external comparator and suggested additional algorithms; expand only for a specific unresolved claim.

Completion criterion: any claimed computational gain is a measured gain for a comparable target, and any cap/prior recommendation is supported by sensitivity results.

### Phase 4: Turn Results into the Revision

- [ ] Add a simulation section organised by the three benefit questions, not a catalogue of plots. Include a scenario/metric table, a performance table with Monte Carlo uncertainty, and a compact uncertainty/location figure.
- [ ] Put the computational comparison, diagnostics, and cap/prior sensitivity in a concise appendix where appropriate. Retain detailed developer experiment records in `JOE revision/simulation/`, outside public package documentation.
- [ ] Recompute any empirical quantities affected by recent sampler, likelihood, or mode-finding fixes before retaining the current numerical claims.
- [ ] Explain unconditional mode-location probabilities separately from probabilities conditional on the number of modes; include the two-components/one-mode example.
- [ ] Strengthen the loan and inflation motivation using a concrete inferential question. Treat multimodality as a distributional feature, not automatic evidence of latent groups or multiple economic equilibria. Decide whether the DNA application belongs in the appendix or should be removed, as the response draft proposes.
- [ ] Update the abstract, introduction, conclusions, and editor/referee replies only with results actually obtained. Clarify prior choices and algorithm details; distinguish the contribution from established shifts, augmentation, and sparse finite mixtures.

Completion criterion: each retained benefit claim has a corresponding numerical result and limitation, and each simulation/efficiency referee comment points to a specific table, figure, or explicit scope decision.

## Immediate Next Step

Continue Phase 2 from [simulation/results/phase2/calibrate/results.csv](simulation/results/phase2/calibrate/results.csv). Address the seven combined calibration failures on the existing datasets before `--run`. The 600 datasets, protocol, paired-summary code, and resumable runner are ready; preparation and calibration are not production evidence of superiority.