# Instructions for Claude Code: rerun the Referee 1 simulations in R

You are helping Paul Labonne, a co-author of "Flexible Negative Binomial Mixtures for Credible Mode Inference in Heterogeneous Count Data from Finance, Economics and Bioinformatics" (Cross, Hoogerheide, Labonne and van Dijk), which is under revision at the *Journal of Econometrics*. Jamie Cross ran the simulations below in R on his machine. Your job is to:
1. rerun them on Paul's machine;
2. check that the results reproduce;
3. help Paul compare them with his own code.

**This folder is self-contained.** It is the unzipped `Pilot_R1.zip` and holds everything you need: this file, `R1_pilot_checks_for_Paul.md`, `code/`, `results/` and `figures/`. Paul has no access to Jamie's files. Do not look for, or assume, any path on Jamie's machine.

Read `R1_pilot_checks_for_Paul.md` first. It explains what each result is for, and which numbers the paper and the replies to the referee quote.

**The simulations (Part A) are the priority.** The real-data runs (Part B) are optional and need the data files on Paul's machine.

## Ground rules

- Work only inside this folder and the output folder you create in it. Do not edit anything outside it, including Paul's own code, unless he asks.
- Do not change the model, the priors, the seeds or the run settings unless Paul asks. The point is to reproduce Jamie's numbers exactly, or to explain why they differ.
- Ask Paul before starting any run that takes more than an hour, and tell him how many cores it will use.
- Record the machine, the R version and the number of free cores, so that timings can be judged. Jamie's timings came from a heavily loaded machine and are not comparable.
- If a script fails, report the error and the fix you propose; do not silently change the method.

## Setup

1. **R.** Use R 4.3 or later; Jamie used R 4.5.3. The packages are base R, `parallel` and, for step A1 only, `quadprog`. Install `quadprog` if it is missing.
2. **Working folder.** Create an output folder, for example `run/` inside this folder, and copy every file from `code/` into it. All scripts read and write relative to one folder `S`:
   - set the environment variable `PILOT_R1_DIR` to that folder's full path, with forward slashes;
   - or run the scripts from inside that folder, since each script uses `S <- Sys.getenv("PILOT_R1_DIR", unset = getwd())`.
3. **Data.** The simulations in Part A need no data files. The default data are typed into `s2_functions.R`. Only Part B, and the real-data tasks inside step A2, need two files that Paul has on his machine:
   - `DATA_d4z4_4_411obs_a1_c2_y3.csv` (DNA);
   - `michigan_15_april_2023.csv` (inflation).

   Ask Paul where they are, or search his machine for these names. Set the environment variable `PILOT_R1_DATA` to their folder, with forward slashes and a trailing slash. Check that the DNA file gives 411 observations and that the stand-in inflation sample (March 2023, responses 0 to 50) gives 518, for example:

   ```
   Rscript -e 'source("snb_s2_sampler.R"); source("s2_functions.R"); sapply(c("default", "dna", "inflation"), function(d) length(real_data(d)))'
   ```

   This should print 4691, 411 and 518.
4. **Smoke test.** Run `Rscript NB_Comparisons.R`. It runs a short example on one simulated data set and prints one line per algorithm. It should finish in about a minute without errors.

Run each long job in the background and log its output to a file, for example `Rscript s2_full.R 6 > full.log 2>&1`. The task scripts skip tasks whose result file already exists, so an interrupted run can be restarted with the same command.

## Part A: simulations (the priority)

Steps A3 and A4 read the results of step A2, so run A2 first.

| Step | Purpose (R1 comment) | Command | Approximate time on Jamie's machine | Output |
|---|---|---|---|---|
| A1 | Limits in Table D.1 (4a) | `Rscript theory_lp_fast.R`, then `Rscript theory_map_fast.R` | minutes | `theory_lp.rds`, printed limits |
| A2 | Shifted against unshifted mixtures, Table D.1 (4a) | `Rscript s2_full.R 6`, then `Rscript s2_full_report.R > full_report.log` | 8 h on 6 workers | `res_full/`, `res_full_J20/` |
| A3 | Shifted NB against shifted Poisson, Table D.2 (4c, 4d) | `Rscript run_4c.R 4`, then `Rscript report_4c.R > report_4c.log` | 16 h on 4 workers | `res_4c/` |
| A4 | Five samplers, Table E.1 and Figure E.1 (4b, 4d) | `Rscript run_nbcomp.R 4`, then `Rscript report_nbcomp.R > nbcomp_report.log`, then `Rscript fig_appE.R` | 30 min on 4 workers | `res_nbcomp/`, `figure_appE.pdf` |
| A5 | Histogram peaks in one NB sample (6) | `Rscript pilot_hist.R > hist.log` | about 1 min | printed table |

**Notes:**
- **Step A2** also fits the three real data sets with 4 chains each. Without `PILOT_R1_DATA`, the DNA and inflation tasks are logged as errors and skipped, while all simulation tasks still run. The default data needed by A4 are typed in, so A4 is not affected.
- **For a spot-check rather than the full runs:** run A1, A4 and A5 in full. For A2 and A3, restrict the `for` loops that build the task list in `s2_full.R` and `run_4c.R` to one or two designs, keeping the seeds, and compare with the reference values below. Tell Paul which designs you chose.

## Part B: real data (optional; needs the data files)

| Step | Purpose (R1 comment) | Command | Approximate time | Output |
|---|---|---|---|---|
| B1 | Real data, both models (3b, 3c, 5, 1a, 6, 7, 9b) | `Rscript pilot_real.R 4`, then `Rscript pilot_report.R > pilot_report.log`, then `Rscript pilot_figs.R` | 45 min on 4 workers | `res_pilot/`, `pilot_fig_*.png` |
| B2 | Table 3, log-predictive likelihood (minor 7) | `Rscript smc_batch.R ours` and `Rscript smc_batch.R sp` (two processes); then rerun `pilot_report.R` | 1–3 h each | `res_smc/`; summary at the end of `pilot_report.log` |

The inflation results use a stand-in sample. Ask Paul for the 480-observation sample used in the paper (see "After the runs").

## Reference values to reproduce

The seeds are fixed in every script, and each task sets its own seed, so results should match Jamie's exactly under R 4.5.3 regardless of the number of workers. Under another R version, small Monte Carlo differences are possible: percentages within a few points, TVs within .005, WAIC within a few points. Report anything larger. Jamie's logs are in `results/`: `report_4c.log` (A3), `nbcomp_report_final.log` (A4) and `pilot_report.log` (B1, B2).

**A1 (Table D.1, limit column):** .040 (shifted second component), .367 (narrow peak at 50), .112 (three close peaks), .373 (heaping), .367 (rounded normals), 0 for the three controls. Appendix D.1 also quotes .275 for NB(μ = 2, r = 10) shifted by 5, from `theory_map_fast.R`.

**A2 (Table D.1).** Percentage of 50 data sets with the correct most probable number of modes, shifted (S) and unshifted (U), and mean TV.

| Design | n = 2,000: S / U correct | S / U TV | n = 500: S / U correct | S / U TV |
|---|---|---|---|---|
| One negative binomial | 100 / 100 | .011 / .011 | 100 / 100 | .021 / .021 |
| 60% zeros and a second peak | 100 / 100 | .016 / .015 | 100 / 100 | .028 / .028 |
| Two components far apart | 100 / 100 | .020 / .018 | 100 / 100 | .040 / .038 |
| Shifted second component | 100 / 100 | .017 / .062 | 98 / 100 | .038 / .074 |
| Narrow peak at 50 | 100 / 100 | .017 / .377 | 100 / 100 | .035 / .387 |
| Three close peaks | 100 / 0 | .025 / .124 | 100 / 0 | .047 / .141 |
| Heaping at multiples of five | 78 / 0 | .031 / .377 | 18 / 0 | .061 / .384 |
| Rounded normal components | 100 / 100 | .047 / .371 | 90 / 100 | .078 / .374 |

Real data inside A2, four chains (WAIC, shifted / unshifted): default 13,974 / 14,007; DNA 3,348 / 3,347; inflation stand-in 2,776 / 3,010.

**A3 (Table D.2).** Correct % (NB / P), TV (NB / P), and the number of data sets out of 50 with lower WAIC for NB.

| Design | n = 2,000 | n = 500 |
|---|---|---|
| One negative binomial | 100/80, .011/.026, 50 | 100/50, .021/.049, 47 |
| 60% zeros and a second peak | 100/86, .016/.021, 43 | 100/74, .028/.040, 45 |
| Two components far apart | 100/0, .020/.040, 50 | 100/0, .040/.075, 48 |
| Shifted second component | 100/92, .017/.029, 47 | 98/94, .038/.071, 48 |
| Narrow peak at 50 | 100/100, .017/.029, 49 | 100/80, .035/.067, 49 |
| Three close peaks | 100/94, .025/.039, 49 | 100/100, .047/.063, 44 |
| Heaping at multiples of five | 78/60, .031/.033, 40 | 18/36, .061/.066, 36 |
| Rounded normal components | 100/100, .047/.044, 15 | 90/98, .078/.071, 8 |
| Three shifted Poisson | 100/100, .026/.021, 4 | 100/100, .051/.041, 1 |
| Two shifted Poisson | 100/100, .019/.016, 6 | 100/100, .037/.031, 10 |

**A4 (Table E.1).** Medians over 24 chains:

| Sampler | IF log-lik | IF modes | IF dispersion | ESS/s log-lik | ESS/s dispersion | CPU s | Chains far from long run |
|---|---|---|---|---|---|---|---|
| ours | 2.0 | 1.3 | 5.8 | 54 | 25 | 15.4 | 1 |
| pr_rw | 9.8 | 1.7 | 138 | 11 | 1.0 | 16.1 | 1 |
| zhou | 7.5 | 2.3 | 31 | 15 | 2.6 | 18.9 | 1 |
| zhou_mh | 5.1 | 1.4 | 33 | 20 | 2.8 | 18.4 | 5 |
| mur_rw | 1.8 | 1.5 | 3.0 | 57 | 41 | 15.4 | 2 |

The inefficiency factors should reproduce. ESS/s and CPU depend on the machine.

**A5 (comment 6).** Two or more peaks, and two or more peaks each holding at least 1% of the sample, with the strict rule in brackets:
- n = 270: 99.2% (97.8%) and 95.7% (93.4%);
- n = 411: 98.8% (97.2%) and 88.9% (86.6%);
- n = 4,691: 91.5% (87.6%) and 5.0% (5.0%).

**B1 (real data):**

| Data, model | Pr(number of modes) | Modes with Pr ≥ .1 | WAIC | Acceptance r̃ / α |
|---|---|---|---|---|
| Default, shifted NB | 2: .94, 3: .06 | 0, 4 | 13,974 | .37 / .54 |
| Default, shifted Poisson | 2: .16, 3: .72, 4: .12 | 0, 4, 29–31 | 13,985 | – / .45 |
| DNA, shifted NB | 1: .90, 2: .10 | 19–23 (top 20) | 3,348 | .28 / .66 |
| DNA, shifted Poisson | 4: .31, 5: .45, 6: .15 | 12–13, 23–24, 38–41 | 3,353 | – / .45 |
| Inflation stand-in, shifted NB | 8: .31, 9: .64 | 0, 3, 5, 10, 15, 20, 25, 30, 50 | 2,779 | .36 / .38 |
| Inflation stand-in, shifted Poisson | 7: .10, 8: .61, 9: .29 | 0, 3, 5, 10, 15, 20, 26–28, 30, 50 | 2,799 | – / .38 |

**B2 (Table 3, log-predictive likelihood):**

| Data | Shifted NB, seeds 1–3 | Shifted Poisson, seeds 1–3 |
|---|---|---|
| Default | −7,016.9, −7,020.0, −7,027.2 | −7,034.2, −7,033.5, −7,038.2 |
| Inflation stand-in | −1,489.7, −1,492.8, −1,496.3 | −1,489.5, −1,487.8, −1,487.9 |
| DNA | −1,683.5, −1,683.3, −1,683.7 | −1,707.6, −1,708.7, −1,708.2 |

On the default data our three seeds spread by about 10, close to the gap of 14 between the models. Once the reference runs reproduce, rerun the default data with more particles, for example `Rscript smc_ml.R default ours 4000 10 1` and the same for `sp`. Report the spread across three seeds.

## After the runs

1. **Write `run/REPRODUCTION_REPORT.md`.** For each step, give the reference value, your value and whether they agree within the tolerances above. Record the R version, the machine, the cores used and the run times. List every deviation, with your diagnosis.
2. **Help Paul check the same quantities in his own code.** This is the main purpose. Use the items in `R1_pilot_checks_for_Paul.md`, simulations first:
   - Tables D.1, D.2 and E.1 on a few designs;
   - acceptance and mixing of the r̃ and α steps;
   - conditional mode locations;
   - the empirical statements;
   - Table 3 and its Monte Carlo error.

   Present each comparison as a table, Paul's code against the R pilot, and flag where they disagree. Do not modify Paul's code unless he asks.
3. **If Paul supplies the 480-observation inflation sample:**
   - add a loader for it to `real_data()` in `s2_functions.R`, keeping the stand-in as an option;
   - rerun B1 and B2 for inflation;
   - report the number and location of modes for both models, so the statements in Section 4.2 of the paper can be checked.
4. **Report back to Paul,** who will pass the numbers to Jamie.

## Files in `code/`

| File | Role |
|---|---|
| `NB_Comparisons.R` | Sampler with a switch between updates: `ours`, `mur_rw`, `pr_rw`, `zhou`, `zhou_mh`, `sp` |
| `snb_s2_sampler.R` | Sampler used for step A2 (shifted and unshifted NB); also supplies `simulate_snb` |
| `s2_functions.R` | Designs G1–G8, true pmfs, `chain_init`, `summarize_chain`, `split_rhat`, `real_data`, `run_task` |
| `theory_lp_fast.R`, `theory_map_fast.R` | A1 |
| `s2_full.R`, `s2_full_report.R` | A2 |
| `run_4c.R`, `report_4c.R` | A3 |
| `run_nbcomp.R`, `report_nbcomp.R`, `fig_appE.R` | A4 |
| `pilot_hist.R` | A5 |
| `pilot_real.R`, `pilot_report.R`, `pilot_figs.R` | B1 |
| `smc_ml.R`, `smc_batch.R` | B2 |
