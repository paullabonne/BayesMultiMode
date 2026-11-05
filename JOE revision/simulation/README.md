# Revision Simulations

This study assesses reconstruction of the overall mixture PMF, its predictive probabilities, and its modes. It does not assess recovery of individual components, labels, parameters, or allocations.

## Run

From the package root, using the normal R installation:

```bash
Rscript "JOE revision/simulation/run-phase1.R" --checks-only
Rscript "JOE revision/simulation/run-phase1.R" --smoke
Rscript "JOE revision/simulation/run-phase1.R"
Rscript "JOE revision/simulation/diagnose-phase1.R" --sparsity
Rscript "JOE revision/simulation/run-phase2.R" --prepare
Rscript "JOE revision/simulation/run-phase2.R" --smoke
Rscript "JOE revision/simulation/run-phase2.R" --calibrate
Rscript "JOE revision/simulation/diagnose-phase1.R"
```

The runner also works from another directory when given its full script path. It resolves package loading and result paths relative to the script, not the working directory. Do not source the runner interactively.

For the Phase 1 runner, `--checks-only` runs correctness and equivalent-posterior feasibility gates without the mixture pilot; `--smoke` is a reduced execution check. Its unflagged command reproduces the historical pilot design, not the intended-prior Phase 2 protocol. Old results are preserved with their original code hashes. Do not overwrite them after a sampler change.

The Phase 2 runner defaults to preparation only. `--prepare` checks the six population PMFs and saves 600 training/test datasets. `--smoke` checks execution with 36 short chains. `--calibrate` runs one dataset per scenario with four chains per method at the candidate production length. `--run` requests the 100-replication comparison only when calibration has passed with the same protocol and code hashes. A failed calibration cannot be bypassed by running `--run`.

The full runner requires the imports declared in [../../DESCRIPTION](../../DESCRIPTION), plus `pkgload` and `testthat`. It does not install dependencies, alter R's PATH, disable certificate checks, or substitute a different sampler when imports are missing.

## Files

- [phase1.R](phase1.R): settings, random count generation using `rpois()`/`rnbinom()`, certified PMF grids, independent truth, mixture-level metrics, and traditional-coordinate NB updates.
- [test-phase1.R](test-phase1.R): independent mode/boundary/plateau, PMF invariance, scoring, finite-grid balance, quadrature, prior transformation, and initialisation checks.
- [run-phase1.R](run-phase1.R): gated execution, package loading, saved seeds/data/chains/failures, and reports.
- [diagnose-phase1.R](diagnose-phase1.R): targeted longer-chain checks, mixture-function trace PDFs, and the pilot assessment. It uses the saved pilot data, not newly generated favourable datasets.
- [phase2.R](phase2.R): fixed scenario parameters, intended-prior settings, independent seeds/data, and paired comparisons that retain failures.
- [run-phase2.R](run-phase2.R): preparation, smoke test, calibration, and gated main comparison.

The traditional NB update uses a random walk on size `r > 1` with density proportional to `1/(1+r)^2`, exactly the size prior induced by uniform `rtilde` on `(0, 0.5)`. It retains the proposed sampler's augmented mean update and full-data mean bound. Its isolated namespace adapter reuses the existing unshifted-NB mixture sampler without changing public APIs. Validate both single-NB samplers against quadrature and run the fixed-complexity adapter before treating the baseline as feasible. This establishes a comparable target, not a speed advantage.

## Gates and Results

Results are saved under `results/phase1/`, or `results/smoke/` for the reduced run. Generated outputs are ignored by Git. `gates.csv` and `status.md` distinguish passed, failed, blocked, and incomplete work. A failed or blocked gate stops the next stage; inspect saved errors rather than dropping datasets, changing `K`, or accepting numerical fallback allocations.

The manifest records settings, the Git commit, and code hashes. Completed fit groups resume only with the same manifest. Archive an earlier result directory before rerunning after code/settings changes. The data seed for scenario index `s` and replication `r` is `20261105 + 100000*s + 1000*r`; chain seeds add `100*m + c` for method index `m` and chain `c`. Training and test observations are independent RNG draws, and all seeds and generated datasets are saved.

The original pilot used `a0 = 10, A0 = 50`, with mean concentration 0.2. That differs from the intended `a0 = 1, A0 = 200`, with mean 0.005, and favours less sparse weights. Phase 2 uses the intended prior. The Poisson rate priors remain `l0 = 2, L0 = 1` in both comparisons to avoid confounding the concentration-prior check.

Constant mixture summaries are reported separately because R-hat/ESS are undefined for constants; an unchanged mode indicator is not proof that unvisited alternatives have zero posterior probability. Varying summaries require R-hat <= 1.01 and bulk ESS >= 400. Continuous summaries also require tail ESS >= 400; undefined tail ESS for discrete indicators is reported, not used as a gate. Group results now separate the mean/zero-mass checks from modal checks, and report between-chain predictive-PMF total variation and mode-count probability MCSE. Warnings, non-finite scores, invalid allocations, or unresolved quadrature disagreement do not authorise production claims.

The runner uses up to four local R workers, bounded by available physical cores, and closes them when the run finishes or fails. Serial and parallel execution were checked to produce identical draws for the same seeds. PMF scoring uses `matrixStats` log-sum-exp routines; `matrixStats` is installed through the declared `posterior` dependency.

All MCMC diagnostics and performance measures concern label-invariant mixture functions. Mode locations use the leftmost point of each adjacent plateau. The grid extends beyond every component's rightmost mode and separately controls omitted tail mass. A genuine zero predictive probability retains its infinite log-score penalty and is not floored.

## Phase 1 Closeout

Phase 1 is closed as a pilot assessment, not as a claim that all samplers have converged. All imports are installed; independent correctness and equivalent-posterior baseline checks pass. The original 480-chain pilot and two longer-chain checks remain archived. Its 29/120 combined pass count must not be interpreted as a diagnosis that entire mixture distributions fail to mix.

The intended prior exposed rounded-zero Dirichlet weights in the old weight-conditioned concentration update. Count samplers now update concentration conditional on allocation counts with weights integrated out, immediately followed by the conditional weight draw. This retains the model's posterior target and matches the manuscript's partition-conditioned concentration formula. Tests check the integrated target against independent beta integration and enforce update ordering. No positive floor is imposed on mixture weights.

The matched follow-up reran both priors with that same corrected sampler: 96 chains, one saved dataset per scenario, four methods, four chains, 3000 iterations / 1500 burn-in. Neither prior caused fitting errors. Both had 4/12 groups passing the mean/zero-mass checks and 2/12 passing the combined checks. Under the intended prior, maximum pairwise between-chain PMF TV ranged from 0.0007 to 0.064. Thus the prior mismatch was real, but it did not explain all diagnostic failures. See [results/phase1-sparsity/sparsity-comparison.csv](results/phase1-sparsity/sparsity-comparison.csv); each saved group summary also contains mode probabilities and their MCSE.

The code/manuscript audit confirms the bounds `mu in (0, max(data))`, `rtilde in (0, 0.5)`, and integer shifts from zero to the full-data maximum, truncated by assigned observations for nonempty groups. Burn-in alignment is explicit in the study. Remaining manuscript work: describe the code's local shift MH update rather than an exact grid draw; distinguish the success-probability convention in the text from the failure-probability variable in code; and correct the likelihood-surface example's `r = 1` versus `r = 3` inconsistency. Detailed prior-sensitivity and speed claims remain future evidence, not established by the audit.

## Phase 2 Status

The six-scenario protocol is implemented with 1000 training observations, 2000 test observations, three methods, and an initial budget of 100 independent repetitions. The seed schedule is separate from Phase 1: `20261106 + 1000000 * scenario_index + 1000 * replication`, followed by `100 * method_index + chain` for fitting. All methods share each generated dataset. No external structural-zero model is added to this initial three-method design; its referee response remains a separate scope decision.

Preparation generated all 600 datasets and the [population PMF figure](results/phase2/prepare/population-pmfs.png). The 36-chain smoke test and 72-chain calibration both completed without fit or scoring errors. Calibration used 12000 iterations / 6000 burn-in: 14/18 groups passed the mean/zero-mass checks, 12/18 passed modal checks, and 11/18 passed the combined gate. These are diagnostic results, not comparisons averaged over independent simulated datasets.

See [results/phase2/calibrate/results.csv](results/phase2/calibrate/results.csv) and [results/phase2/calibrate/status.md](results/phase2/calibrate/status.md). Phase 2 is underway, but the full 7200-chain comparison has not started: remaining calibration failures are retained, and the production guard remains active. Paired estimates require multiple independent datasets and retain non-finite/failed pairs explicitly; unsuccessful fits are never silently removed.