# Referee 1: pilot results and checks for the main code

*Draft, 9 November 2026. Prepared by Jamie (with Claude) for Paul.*

## Purpose

The replies to Referee 1 now have text for every comment. Several replies and manuscript passages still depend on numbers from the main code (in your private repository). To make the remaining work concrete, we ran a pilot of every outstanding R1 computation in R. This note gives, item by item:

- what the paper or the reply currently says;
- what the R pilot gives;
- what to check or produce in the main code.

The pilot is a reference point, not the final result. The paper's numbers should come from the main code. Where the main code and the pilot disagree, please tell us which is right and why.

**The simulations are the priority:** Section 8 below (Tables D.1, D.2 and E.1, for R1 comment 4), and the histogram simulation in Section 7. The real-data items in Sections 1–6 come second.

This folder (the unzipped `Pilot_R1.zip`) holds everything: this note, the R code, Jamie's result logs and the pilot figures. `CLAUDE_INSTRUCTIONS_R1_pilot.md` can be given to Claude Code to rerun the R simulations on your machine, check them against our numbers, and compare them with your code.

## How the pilot was run

- **Code.** `NB_Comparisons.R`: one sampler in base R, with a switch for the update of the component parameters. `"ours"` is the sampler of Section 3.1; `"sp"` is the shifted Poisson mixture, our model with r̃_j = 0. Everything else is shared: allocation, joint κ step, empty components drawn from the prior, α step, weights.
- **Settings, as in the paper.** J = 10. Sparse prior α ~ Gamma(1, rate 200). μ_j ~ U(0, M], κ_j uniform on {0, …, M}, r̃_j ~ U(0, ½), with M = max(y).
- **Run length and chains.** 50,000 burn-in and 50,000 kept iterations; 4 chains per model and data set.
- **Starting values.** Chain 1 starts at r = 10 and κ = 0. Chains 2–4 start at r = 2 with random shifts.
- **r̃ random walk.** It starts at sd 0.1 and is adapted in batches of 50 iterations during burn-in towards acceptance 0.3–0.4, then fixed.
- **α step.** A random walk on log α with step s = 1, not adapted.
- **Inflation data.** We do not have the paper's 480-observation sample, so the pilot uses a stand-in: March 2023, responses 0 to 50, n = 518. The inflation results below are therefore illustrative only.
- **Timings.** The machine was heavily loaded by other jobs, so CPU times in this note are not comparable. Please time on a quiet machine.

The scripts are listed at the end.

---

## 1. Comment 3b: acceptance and mixing of the r̃ step

**Reply says.** The candidate is Gaussian, and candidates outside (0, ½) are rejected. The sd of each r̃_j is adapted in batches during burn-in (Roberts and Rosenthal, 2009), then fixed, targeting acceptance 0.3–0.4. Red note: *"Paul to implement this adaptation, and to add the acceptance rates and mixing, measured on label-invariant quantities."*

**Pilot: real data.** Pooled over non-empty components and 4 chains.

| Data | Acceptance of r̃ | Adapted sd (median; range) | IF of Σπ_j log r_j | IF of log-likelihood |
|---|---|---|---|---|
| Default | 0.37 | 0.34; 0.03–0.49 | 8.6 | 11.2 |
| Inflation (stand-in) | 0.36 | 0.49; 0.06–0.56 | 165 | 324 |
| DNA | 0.28 | 0.28; 0.04–0.42 | 69 | 650 |

IF is the inefficiency factor per stored draw (every 5th of 50,000), median over chains.

**Pilot: simulations.** In the five-sampler comparison (Appendix E), the median acceptance of r̃ over components is 0.40–0.50 by design (overall range 0.33–0.62). That is somewhat above the 0.3–0.4 target. The adaptation is per component label, and components empty and refill after burn-in, so the fixed sd need not suit the component that later occupies a label.

**Please check in the main code:**
1. Is the adaptation implemented as the reply describes (batches during burn-in, target 0.3–0.4, then fixed)? If not, either implement it or change the reply to what the code does.
2. Report acceptance and the IF of label-invariant quantities for the three data sets: log-likelihood, number of modes, Σπ_j log r_j over non-empty components, and α. A small table in Appendix C would answer the referee.
3. The mixing on the DNA and inflation data is slow. Log-likelihood IFs of about 300–650 mean an effective sample of only about 15–30 from 10,000 stored draws per chain. Please check whether the main code shows the same. If it does, we should use longer runs for the applications, or say so in the paper.

## 2. Comment 3c: the α step

**Reply says.** a_α = 1 and b_α = 200, with a random walk on log α (Section 3.1). Red note: *"Paul to add s, the acceptance rate and mixing."*

**Pilot (s = 1).**

| Data | Shifted NB: acceptance | Shifted NB: IF of α | Shifted Poisson: acceptance | Shifted Poisson: IF of α |
|---|---|---|---|---|
| Default | 0.54 | 1.8 | 0.45 | 31 |
| Inflation (stand-in) | 0.38 | 14 | 0.38 | 8.6 |
| DNA | 0.66 | 68 | 0.45 | 3.9 |

**Please check:** the value of s in the main code, its acceptance rate and the IF of α on the three data sets, and add s to Section 3.1.

## 3. Comment 5: mode locations conditional on the number of modes

**Reply says.** Figures 5–7 now show Pr(y is a mode | L = l) whenever the posterior of L puts mass on more than one value. Red note: *"Paul to add the conditional location plots to Figures 5–7."*

**Pilot.** The figures `pilot_fig_default.png`, `pilot_fig_inflation.png` and `pilot_fig_dna.png` show one way to do this. The right-hand panels give the marginal Pr(y is a mode) as grey bars and the conditional probabilities for the two most probable values of L as points. Examples:
- **Default, shifted NB.** Given L = 2 (Pr .94), the modes are at 0 and 4 with probability 1. Given L = 3 (Pr .06), there is an added mode around 29 (.21).
- **DNA, shifted Poisson.** Given L = 5 (Pr .45), the modes are near 12–13, 23–24 and 39–40, plus low-probability modes further out.

**Please produce:** the conditional panels in the main code's Figures 5–7, in whatever form you prefer.

## 4. Minor comments 2 and 6: figure formatting

**Requests:**
- larger axis fonts and thicker lines in Figures 2–4;
- no black grid in Figure 4;
- in the "Number of modes" plots, at least three values on the horizontal axis and the same grid for both models;
- confirm what the red lines are and how many draws are shown.

**Pilot.** The pilot figures use one common grid of modes for both models, larger fonts, thicker lines and no grid. The red lines are the probability mass functions of 400 posterior draws (100 per chain), with the posterior mean in black.

**Please:**
1. Update the figures.
2. Tell us what the red lines are in the main code's figures and how many draws they show. The reply currently says "fitted probability mass functions for posterior draws".
3. Write the two short replies for minor comments 2 and 6 once the figures are done.

## 5. Empirical statements to recheck

Comments 1a, 6, 7 and 9b, Sections 4.1–4.3 and Section 5 quote the results below. The pilot agrees on the main points but not on several details.

| Statement in the paper or replies | Pilot | Action |
|---|---|---|
| Default, shifted NB: "100% posterior probability for two modes" (Section 4.1) | Two modes .94, three modes .06 (the third mode is around 29) | Use the main code's value. "100%" is probably too strong. |
| Default, shifted NB: modes at 0 and 4 | 0 and 4, each with probability 1.00 | Agrees |
| Default, shifted Poisson: "puts quite a bit of mass … in the far right tail around 10–20"; "more than two modes (albeit with low probability)" | Three modes .72; tail modes near 29–31 (and 17–18 given L = 4) | Recheck the location of the tail mode and the probability of more than two modes. "Low probability" does not fit .84. |
| DNA, shifted NB: "posterior probability of unimodality is 100%"; "mode between about 17 and 21, highest at 19" | Unimodal .90; mode probability spread over 19–23, highest at 20 (.21) | Recheck both. "100%" and the location may change. |
| DNA, shifted Poisson: "four modes, at about 12, 22, 38 and 62" | Five modes .45, four .31; modes near 12–13, 23–24, 38–41, and low-probability modes near 55–60 | Recheck. "Four modes" may need to become "four or five". |
| WAIC (4a reply): default 13,974 against 14,007 (unshifted); DNA 3,348 against 3,347; inflation 2,776 against 3,010 | Shifted NB against shifted Poisson: default 13,974 against 13,985; DNA 3,348 against 3,353; inflation (stand-in) 2,779 against 2,799 | Confirm the shifted NB values in the main code (the unshifted values are from the earlier R runs) |
| Inflation: unimodal 2%; modes at 3 and 5; shifted Poisson modes near 3 and 20–24 | Stand-in sample: 8–9 modes for both models (0, 3, 5, 10, 15, 20, 25, 30, 50) | Cannot be checked without the 480-observation sample. Please send it, or tell us how it is built from the Michigan file (month, range, handling of zero and negative responses). |

**Convergence warnings from the pilot:**
- **Shifted Poisson:** split R̂ of the log-likelihood is 1.53 on the default data and 1.47 on the stand-in inflation data. The chains disagree in level, although all four agree on the most probable number of modes.
- **Shifted NB:** R̂ is at most 1.05 on all three data sets, but mixing on DNA and inflation is slow (Section 1 above).

## 6. Minor comment 7: Table 3, the log-predictive likelihood

**Paper says.** Table 3 notes: "computed from the sequence of one-step-ahead predictive likelihoods of all observations (Geweke, 2005)". Red note: the current values (default −6,984 against −7,015; inflation −1,373 against −1,374; DNA −1,674 against −1,773) are the posterior mean of the in-sample log-likelihood and must be recomputed.

**Pilot method.** The sum of one-step-ahead log predictive likelihoods equals the log marginal likelihood. We estimate it by sequential Monte Carlo that adds one observation at a time (Chopin, 2002), so each incremental weight is exactly p(y_t | y_1, …, y_{t−1}). Details:
- 1,000 particles drawn from the prior.
- Resampling when the ESS falls below 500.
- 5 sweeps of the sampler per move on the data seen so far.
- Three seeds, each with a different random order of the observations, to show the Monte Carlo error.

**Pilot values.** Three seeds each; the mean is in brackets.

| Data | Shifted NB | Shifted Poisson | Difference in means | Current Table 3 (SNB / SP) |
|---|---|---|---|---|
| Default | −7,016.9, −7,020.0, −7,027.2 (−7,021.4) | −7,034.2, −7,033.5, −7,038.2 (−7,035.3) | +14 for SNB | −6,984 / −7,015 |
| Inflation (stand-in) | −1,489.7, −1,492.8, −1,496.3 (−1,492.9) | −1,489.5, −1,487.8, −1,487.9 (−1,488.4) | −4 for SNB | −1,373 / −1,374 |
| DNA | −1,683.5, −1,683.3, −1,683.7 (−1,683.5) | −1,707.6, −1,708.7, −1,708.2 (−1,708.2) | +25 for SNB | −1,674 / −1,773 |

**What the pilot suggests:**
- **DNA:** our model is clearly ahead, by about 25 log points, with small Monte Carlo error (seeds spread by at most 1.1). The paper's gap of 99 comes from comparing in-sample averages.
- **Default:** our model is ahead by about 14 points, and every seed of ours lies above every seed of shifted Poisson. But our seeds spread by 10, so the size of the gap is uncertain. With n = 4,691, 1,000 particles and 5 sweeps are not enough. Please use more particles or sweeps, or another estimator, and report the Monte Carlo error.
- **Inflation (stand-in):** shifted Poisson is slightly ahead, by about 4 points, within the spread of our seeds (6.7). The current table also does not separate the models. The paper's sample is needed before anything can be said.
- **The text will need to change.** The Introduction says "a good or better fit … in all three datasets", and Section 4.3 says the log-predictive likelihood favours our model for DNA. Both should be rechecked against the main code's values.

**Please consider:**
1. **How the main code computes Table 3.** If it uses another estimator of the marginal likelihood (for example bridge sampling, or Chib's method adapted to label switching and the discrete shifts), please compare it with the pilot on the DNA data.
2. **Dependence on M.** Unlike WAIC, the marginal likelihood depends on the prior. The bound M = max(y) sets the uniform ranges of μ_j and κ_j, and our model has the extra prior on r̃_j. Please use the same M for both models, and rerun with a larger M as a check. This links to the manuscript's red note *"rerun with different large values of M"*.

## 7. Comment 6: histogram peaks in a single negative binomial sample

**Reply says.** For NB(μ = 5, r = 5), the histogram shows at least two peaks in about 99% of samples of size 270, and in about 95% when values holding under 1% of the sample are ignored. With n = 4,691 such peaks occur in about 5% of samples, while a further peak in the tail occurs in about 92%. Red note: *"Paul: reproduce this simulation … The DNA sample may be 411 alleles rather than 270, in which case use n = 411."*

**Pilot.** 10,000 samples per size. A peak is a value, or a run of equal counts, higher than both neighbours. The "strict" rule counts only single values.

| n | Two or more peaks | … each holding ≥ 1% of the sample |
|---|---|---|
| 270 | 99.2% (strict 97.8%) | 95.7% (strict 93.4%) |
| 411 | 98.8% (strict 97.2%) | 88.9% (strict 86.6%) |
| 4,691 | 91.5% (strict 87.6%) | 5.0% (strict 5.0%) |

The pilot reproduces the reply's numbers.

**Please confirm** the DNA sample size (270 or 411). If it is 411, the reply should say about 89% rather than about 95%.

## 8. Comment 4 (a–d): numbers to confirm in the main code

These are already in the replies and in Appendices D and E, all from the R simulations.

- **Table D.1 (4a):** shifted against unshifted mixtures, eight designs, 50 data sets at n = 500 and 2,000, and the limits from the linear program.
- **Table D.2 and the efficiency paragraph in D.4 (4c, 4d):** shifted NB against shifted Poisson, ten designs.
- **Table E.1 and Figure E.1 (4b, 4d):** five samplers of the component update.

This is the priority. You don't need to rerun everything in your own code: a spot-check of a few designs with your sampler, against the R values in the tables, would be enough to confirm that the two implementations agree. The Claude Code instructions list the reference values design by design.

## 9. Other checks raised in the R1 replies

- **Minor comment 5.** Confirm that the main code draws κ_j given m_j = κ_j + μ_j exactly, over the whole grid {0, …, min(min y in the component, M)} with (0, M] for m_j − κ_j. The R code does.
- **Comment 1b.** A literature check, before submission, for earlier uses of r̃ and for negative binomial mixtures with unknown shifts (Jamie and Lennart).
- **Timings.** Please time the shifted NB and shifted Poisson mixtures, and the five samplers of Appendix E, on a quiet machine.

## Files

Everything is in this folder:

- **`code/`**
  - `NB_Comparisons.R`: the sampler, with the switch between updates.
  - `snb_s2_sampler.R`, `s2_functions.R`: the simulation designs, starting values, summaries and data loading. They are sourced by the scripts below.
  - `pilot_real.R`, `pilot_report.R`, `pilot_figs.R`: the real-data chains, the report and the figures.
  - `smc_ml.R`, `smc_batch.R`: the Table 3 pilot.
  - `pilot_hist.R`: the comment 6 simulation.
  - `run_4c.R`, `report_4c.R`, `run_nbcomp.R`, `report_nbcomp.R`, `fig_appE.R`: the Appendix D.4 and E results.
  - `s2_full.R`, `s2_full_report.R`: the 4a run (Table D.1), which the D.4 and E reports also read.
  - `theory_lp_fast.R`, `theory_map_fast.R`: the limits in Table D.1 and the example in Appendix D.1.
- **`results/`:** `pilot_report.log` (real data and Table 3), `report_4c.log` (Appendix D.4), `nbcomp_report_final.log` (Appendix E).
- **`figures/`:** `pilot_fig_default.png`, `pilot_fig_inflation.png`, `pilot_fig_dna.png`.

**`CLAUDE_INSTRUCTIONS_R1_pilot.md`** is a set of instructions you can give to Claude Code to rerun everything in R and check it against our numbers. It covers setup, the order of the runs, run times, reference values with tolerances, and the report to write. The scripts read their working folder from the environment variable `PILOT_R1_DIR`, or use the folder they are run from, and the data folder from `PILOT_R1_DATA`. The `.rds` result files are not copied, because they are large; we can share them if needed.
