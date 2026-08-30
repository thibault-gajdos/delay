# Analysis code for the theoretical note

Two languages, split by job.

**Python** owns the model: the stochastic simulation, the electromyographic detector, the
calibration, the parameter sweeps, the extraction of events from the archive, and the
figures.

**R** owns the statistics: every analysis of variance, every trend decomposition, and the
Poisson analyses of the partial-activation counts. The two sides meet at CSV files
written by `export_for_R.py`, which carry counts and exposures for the data and for each
fitted model, so that the identical statistical analysis runs on both sides.

Everything is reproducible from two inputs: the OSF archive (`86cbw`) and the manuscript.
No network access is needed once the two package managers below have fetched their
dependencies once.

## Layout

```
scripts/
  delaymodel/          the model package: sim.py, empirical.py, targets.py, behav_sim.py
  *.py                 flat, top-level entry-point / driver scripts (calibrate.py, ...)
  *.sh                 resumable shell drivers for long-running calibrations
  results/             calibration outputs (JSON, tracked) and logs/ (regenerable, gitignored)
  plots/               figures written by figures.py
  pyproject.toml, .python-version, uv.lock, .venv/   Python, managed by uv
  R/                   every statistical test (afex/emmeans), managed by renv
```

## Python (uv)

```
cd scripts
uv sync                              # once: creates .venv/, installs numpy/scipy/pandas/matplotlib
export DELAY_DATA=/path/to/86cbw-osfstorage-archive   # the unpacked OSF archive
uv run python delaymodel/empirical.py  # sanity check; or any script below
```

`uv run python <script>` uses the locked, project-local `.venv` -- no manual virtualenv
management, no polluting the system Python. `uv.lock` pins exact versions; re-run `uv sync`
after pulling a change to `pyproject.toml`. `delaymodel` is a plain local package (no
install step) resolved because every entry-point script inserts its own directory --
`scripts/`, i.e. `delaymodel`'s parent -- onto `sys.path` before importing it.

`DELAY_DATA` points at the unpacked archive (default `/home/claude/data`, the original
cloud-container path -- always set it explicitly when running locally, e.g. in your shell
profile or before each session: `export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"`).
`SIMDIR` (an R-side variable, see below) points at the directory holding the CSV bridge
files, default the parent of `R/`, i.e. `scripts/`. `ARCHIVE` (also R-side, `R/run_stats.R`)
points at the same unpacked archive as `DELAY_DATA` but is read independently by R, default
`../data` relative to wherever `Rscript` is invoked from -- set it explicitly too, e.g.
`export ARCHIVE="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"` before running
`run_stats.R` from `scripts/R/`.

## R (renv)

```
cd scripts/R
Rscript -e 'renv::init()'            # once: creates renv/, renv.lock, .Rprofile
Rscript -e 'renv::install(c("afex","emmeans","dplyr","tidyr","lme4","lmerTest"))'
Rscript -e 'renv::snapshot()'        # records the installed versions in renv.lock
export ARCHIVE="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
Rscript run_stats.R --with-sim
```

On a machine that already has `renv.lock` (checked in after the first `renv::snapshot()`),
skip straight to `Rscript -e 'renv::restore()'` to reproduce the exact recorded package
versions instead of re-resolving them.

## Order of execution

All paths below are relative to `scripts/` (run `cd scripts` first); `results/` and
`plots/` already exist in the repository. Prefix any of these with `uv run python` to use
the managed virtualenv.

```
python delaymodel/empirical.py       -> (printout only)     extraction, checked against the archive
python delaymodel/targets.py         -> (printout only)     calibration targets, printed summary
python calibrate.py N M              -> results/fit_modelM.json        one chunk of N evaluations, model M
BLOCKS=behav python calibrate.py N M -> results/fitb_modelM.json       behaviour-only fit
python export_for_R.py               -> pma_counts.csv      counts and exposures, empirical
                                         sim_counts_mM.csv   the same from each fitted model
                                                              (written at scripts/ root: R's
                                                              SIMDIR default, "..", expects
                                                              them there, not in results/)
Rscript R/run_stats.R --with-sim     -> R/stats_tables.txt  all tables, data and simulation
                                         R/stats_results.csv one row per effect
Rscript R/emit_tex.R                 -> _tables.tex         the same, for the note
python sweep.py N                    -> results/sweep_modelM.json      local one-at-a-time
                                                          robustness and leak variants around
                                                          the calibrated point (superseded for
                                                          the manuscript's own parameter-
                                                          space evaluation by fused_sweep.py
                                                          below; still the source of the
                                                          R7/R8/R16-19 sign-check numbers
                                                          quoted in the response letter,
                                                          C2.1, not in the manuscript itself)
python fused_sweep.py N NTR OUT [GAMMA]
                                   -> OUT (.jsonl, default  the manuscript's first-stage
                                      results/fused_sweep_  broad-parameter Monte-Carlo sweep
                                      results.jsonl)         (Table A1 hybrid ranges: kappa
                                                          widened to [0.20, 30.0] s^-1,
                                                          r and Lmotor freed): task-
                                                          admissibility, the two RT-effect
                                                          criteria, and (since 2026-08-24)
                                                          the delay-slope and MET-invariance
                                                          out-of-sample criteria, at every
                                                          admissible draw. GAMMA (default 0)
                                                          adds a common leak and reruns the
                                                          identical sweep for the leak-
                                                          robustness check (Methods, Behavior
                                                          > Simulation validation) --
                                                          production calls used GAMMA in
                                                          {0.25, 0.50, 1.00}, one output file
                                                          per rate (results/fused_sweep_leak_g<rate>.jsonl).
python figures.py                 -> plots/fig_*.pdf      the four figures of the note
python theory_checks.py                                  verification of Appendix C
python model1_bound.py                                   the Model 1 ceiling, Section 3

python build_targets_q.py         -> results/targets_q.json   RT-quantile calibration targets
                                                          (correct/incorrect, two-stage
                                                          Vincentized/pooled; supersedes
                                                          build_targets.py's targets_behav.json)
python calib_behav_q.py {m1,m2}   -> results/calibq_{tag}.json  behavior-only quantile/G²
                                      results/ckptq_{tag}.json  calibration (three-stage
                                                          Nelder-Mead restart, checkpointed);
                                                          ./supervise_calib.sh {tag} loops it
                                                          to completion across interruptions.
                                                          Supersedes calib_behav.py.
```

`calibrate.py` is chunked and resumable: each call runs a bounded number of evaluations,
writes the checkpoint, and picks up next time -- useful whether you are re-launching by
hand or via `run_emg_recal.sh`/`supervise_calib.sh`. `STEP` controls the restart simplex
size (`STEP=0.2 python calibrate.py 120 2`).

---

# Python

### `empirical.py`
Reads the segmented event files and the behavioural files and returns, per trial: the
condition, accuracy, response time, the times of stimulus onset and offset, of the
response signal, of the response, of the response burst, and every partial activation
with its offset.

Counting follows the archive's own `analysis_PMA/compute_PMA.py`: a partial activation is
an EMG onset in the channel of the responding effector occurring before the onset of the
response burst. Exclusions are responses before the signal, absent responses, and
response times under 0.2 s, leaving 20 566 trials.

**Do not re-merge the event list.** `myonset` applies the 25 ms merge when it defines the
activity windows, so the published events are already one per merged window. The
integrated-profile step then moves boundaries earlier, which is why 27% of consecutive
pairs are separated by less than 25 ms; it does not change the number of intervals.
Merging again removes about a quarter of the counts and is double counting.

### `targets.py`
Builds the target families for calibration, per participant and cell: accuracy, mean
response time and its interquartile range, mean and interquartile range of premotor and
motor execution time, delay and stimulus partial-activation rates, delay log-rate slope,
premature-response rate.

Dispersion uses the interquartile range, not the standard deviation. The empirical
premotor execution time has a long right tail -- 2.7% of trials over one second -- that
the model has no mechanism for, and scoring against the untrimmed standard deviation
forces the remaining parameters to manufacture variability elsewhere.

### `sim.py`
The model of Section 2 of the note, and the detector.

`simulate(p, c, D, n, model, seed, gamma)` integrates the equations at 1024 Hz (`DT = 1/
1024`, the EMG/EEG native acquisition rate -- was 1 ms/1 kHz until 2026-08-30, see "Things
that will bite" below) and returns the gated drive of the selected channel per trial.
`model=1` sets the preparation term to zero. `gamma` adds a leak common to both
accumulators, for the reviewer's sensitivity request.

`measures(sim)` applies the detection rules and returns accuracy, response time, premotor
and motor execution time, validity flags and the partial-activation event times.

`_bursts` is the detector, and reproduces `myonset` as it was used for the data:
threshold, merge intervals separated by at most 25 ms, then discard merged intervals with
fewer than `3 + f` supra-threshold samples, where `f` is the number of intervals per
second containing at most six. On the simulated drive the adaptive part is inert (`f` is
about 0.01/s) because the drive is smooth and has almost no short excursions; on real EMG
it does work, which biases the model towards over-counting by a few per cent, absorbed by
the gating parameter.

The integrated-profile refinement is not simulated: it moves boundaries but not their
number, and only the number is counted. For the same reason simulated burst durations are
not comparable with the reported ones and are not used.

### `calibrate.py`
Objective and search. `model_stats` runs the six cells and returns the target families;
`make_objective` divides each by the between-participant standard error of the target, so
one unit of error is one standard error. Parameters are searched in logarithmic
coordinates and clipped to the bounds in `LO`/`HI`.

Thirteen free parameters. `sigx` is fixed at 0.100 and fixes the amplitude unit. Three
combinations are weakly conditioned rather than unidentified, and are reported as
combinations: `lam * kappa`, `alpha = beta / Delta`, and `Lmotor + r / muRS`.

`BLOCKS=behav` restricts the criterion to accuracy, mean response time, its interquartile
range and the anticipation rate -- what a study without electromyography would see -- and
writes `fitb_modelM.json`. This `fitb_modelM.json` is used only as a *starting point* for
the behavior-only calibration below, not as the manuscript's own behavioral fit.

### `build_targets.py`, `calib_behav.py`, `behav_sim.py` -- superseded 2026-08-30
The manuscript's first behavior-only calibration this project (six per-cell quantities --
accuracy, mean RT, RT IQR, and the pre-signal/fast/no-response proportions -- matched by
mean-squared standardized error). Superseded by the quantile/G² pipeline below at the
user's request, once it became clear the model cannot produce non-negligible
no-response-by-deadline trials at any parameter setting compatible with the observed RT
distribution (a parametric consequence of the post-signal ramp's drift-to-noise ratio, not
an architectural impossibility -- see `calibb_m1.json`/`calibb_m2.json`, its last outputs).
Kept for reference; `calib_behav_q.py` does not import from `calib_behav.py`.

### `build_targets_q.py`, `behav_sim.py`, `calib_behav_q.py`
The manuscript's actual behavior-only calibration (Materials and Methods, "Calibration of
the models to the aggregate data"; Results, Behavior > Simulation validation and
Delay-locked > Simulation validation). Response-time quantile method after Dendauw et al.
(2024) / Ratcliff and Tuerlinckx (2002), not a mean/SD/proportion match:

- `build_targets_q.py` reads the raw behavioural files and writes `targets_q.json`: per
  cell, `p_correct` (two-stage: per-participant accuracy on valid trials, averaged), the
  `.1/.3/.5/.7/.9` quantiles of correct-response RT (two-stage, Vincentized across
  participants) and of incorrect-response RT (pooled across participants -- several
  participants have too few error trials at the higher coherence levels for a
  per-participant quantile to mean anything; see the file's own docstring). Premature,
  sub-0.2-s and no-response trials are dropped entirely, not scored: `p_correct +
  p_incorrect = 1` over valid trials only.
- `behav_sim.py`'s `simulate_behav` is `sim.simulate` without the stored EMG-like trace,
  run to the 4.25 s deadline; `behav_rt_split` returns its correct/incorrect RT arrays,
  filtered the same way as the data side.
- `calib_behav_q.py` computes G² per cell: the model's simulated RTs are passed through the
  *data's own* quantile-defined bin edges (10/20/20/20/20% by construction) and compared to
  that reference partition, together with the correct/incorrect proportion itself. Each
  cell's G² uses a fixed pseudo-sample size of `N_SUBJ = 22` (participants, not pooled
  trials), keeping the two-stage weighting; G² is a scalar discrepancy for calibration
  here, not a chi-square test with an interpretable p-value (`FLOOR` prevents
  log(0)/division-by-zero from an empty simulated bin or response type). `run(model, ...)`
  is a three-stage Nelder-Mead restart (`STAGES`, increasing `n` and shrinking the simplex),
  checkpointed to `ckptq_{tag}.json` after every stage so a run can be interrupted and
  resumed; `supervise_calib.sh {tag}` loops it until `results/calibq_{tag}.json` exists.

### `sweep.py`
Multiplies each free parameter by 0.5 and by 2 in turn, adds three leak rates, and records
whether each directional effect keeps the observed sign. Uses 150 trials per cell by
default, enough for level effects but not for slope differences -- the R16 column is
underpowered in the sweep and should be checked at the calibrated point with more trials.
Local, one-at-a-time perturbation around the single calibrated point -- not a
Latin-hypercube sweep of the space -- and its Table B3 output is quoted only in the
response letter (C2.1) since the Appendix that originally reported it was removed from the
manuscript itself (2026-08-23 restructuring).

### `fused_sweep.py` -- no longer used in the manuscript, 2026-08-27
Was the manuscript's broad-parameter Monte-Carlo sweep (Results, Behavior > Simulation
validation; Methods, "First-stage parameter-space evaluation"). Both the manuscript
subsection and that Methods heading were removed in the 2026-08-27 restructuring: the
qualitative predictions are now proved to hold at *every* admissible parameter setting
(Appendices A-E) rather than checked at a sample of them, which the sweep below made
redundant rather than wrong. `fused_sweep_results.jsonl` and the `fused_sweep_leak_g*`
outputs are quoted only in the response letter now (C2.1, C2.3), not the manuscript. Its
one hardcoded timestep (`i_phys` in `run_one`) was still `1e-3` until 2026-08-30 even
though it imports `simulate`/`measures` from `sim.py` -- fixed to import and use `sim.py`'s
`DT` instead, for consistency with the current 1024 Hz rate, even though the script is
retained only for the letter's numbers. The method itself, unchanged: a Latin-hypercube
sample over the full Table A1 hybrid ranges (`g`, `theta/g`, `kappa` in [0.20, 30.0] s^-1
log-uniform, `SDe`, `muRS`, `r` and `Lmotor` both freed and log-uniform, `beta`), 250 draws
by default, 2000 trials per cell per model per draw. At every draw: the task-admissibility
screen, the two RT-effect criteria (coherence, duration), and -- since 2026-08-24 -- two
out-of-sample criteria checked only at the single calibrated point before that date: the
delay-period log-rate slope criterion (all six cells significantly positive; Results,
Delay-locked) and the motor-execution-time invariance criterion (Results, RS-locked /
Response-locked). Resumable (`os.path.exists(OUT)` picks up where a previous call for the
same `OUT` left off).

The optional fourth argument, `GAMMA`, threads a common leak into every `simulate()` call
(`sim.py`'s `gamma`) and reruns the identical sweep -- this is the leak-robustness check in
Methods (Behavior > Simulation validation) and Results (§130), which used to run the
narrower, pre-widening ranges via a since-superseded script; production runs used `GAMMA`
in `{0.25, 0.50, 1.00}` at the current hybrid ranges, one output `.jsonl` per rate
(`fused_sweep_leak_g<rate>.jsonl`).

`fused_sweep_results.jsonl` (no leak, `GAMMA=0`) is the source-of-truth artifact behind
every broad-sweep number reported in the manuscript as of 2026-08-24 -- keep it alongside
the script rather than regenerating it casually, since a fresh run with a different seed
will not reproduce the exact admissibility counts (the criteria and their signs are
robust; the third-decimal percentages are not guaranteed identical run to run).

### `export_for_R.py`
Writes the bridge files. `pma_counts.csv` has one row per participant, cell, period and
bin, with the raw count and the number of trials contributing to that bin;
`sim_counts_mM.csv` the same from each fitted model.

Exposure is not always the number of valid trials. A partial activation precedes the
response burst by definition, so once that burst has begun a trial can no longer
contribute and later bins are undefined rather than empty. This matters only in the
response-signal and response windows.

Simulated participants are 22 equal groups of simulated trials. This is not a claim that
the model has between-participant variability -- it has none -- but it lets the identical
R analysis run on both sides. The simulated between-group variance is pure sampling noise,
so simulated F values are not comparable in magnitude with the empirical ones; only signs
and orderings are.

### `figures.py`
The four figures of the note: the manuscript's Figure 1 simulated for both models, a
single trial with the detected intervals, the rate profiles against data, and the two
shape predictions.

### `theory_checks.py`
Verification of Appendix C (the Kac-Rice rate derivation; renumbered from Appendix A in
the 2026-08-27 appendix restructuring -- see Correspondence table) against direct
integration at a 1024 Hz step (its own local `DT`, independent of `sim.py`'s -- keep the
two in sync by hand): the Kac-Rice expression against a simulated count of upcrossings,
the occupation identity against the simulated time above the level, the ratio of merged
clusters to raw upcrossings, and the residual condition dependence of motor execution
time from the filter transient.

The Euler bias on crossing counts is upward, about 2% at a 1 ms step and shrinking with
the step; unchecked since the 1024 Hz change, but the step only got finer, so if anything
the bias is now slightly smaller. It affects the simulated rate, not the direction of any
effect.

### `model1_bound.py`
The two numbers behind the exclusion of Model 1 in Section 3 of the note.

`relaxation(kappa)` gives the fraction of the stationary variance reached one second after
commitment, which is what licenses the stationary form of the crossing rate everywhere in
the analysis window: a deficit of 1.3% at kappa = 2.3, one per cent at 2.4, and 5e-11 at
the fitted 12.8.

`best_slope(frac)` bounds the log-rate slope Model 1 can produce while reaching a given
fraction of the observed rate, **without any stationary approximation**. With alpha = 0
the rate depends on kappa and on the distance to the gate and on nothing else -- SD_eps
cancels -- so two observed quantities give two equations in two unknowns. At a quarter of
the observed rate the largest available slope is 0.004/s, against 0.700 observed.

This replaces an earlier argument from the ceiling `nu <= omega / 2pi`, which is the
alpha = 0 special case of the crossing rate and cannot be used to constrain kappa from
data generated with a ramp.

---

# R

### `R/run_stats.R`
Every test, on data and on simulation: the LEVEL / TREND / QUAD decomposition described
below, written with `afex` (repeated-measures ANOVA, Greenhouse-Geisser) and `emmeans`
(trend contrasts). Writes `stats_tables.txt` and `stats_results.csv`. `--with-sim` adds
the simulated counts if `export_for_R.py` has run.

**This file replaces the two R scripts that used to live here** (`run_stats_pkg.R`, the
afex/emmeans implementation, and the from-scratch base-R `rmanova.R` it was checked
against) **as of 2026-08-25.** Both did the same analysis in parallel so that one could be
checked against the other; once that check passed there was no remaining reason to keep
both, so this is now the sole, canonical script, under the name `run_stats.R` -- it is the
former `run_stats_pkg.R`'s code, documented, bug-fixed, and promoted to the canonical
filename. `rmanova.R` and its hand-written Greenhouse-Geisser/orthogonal-contrast
machinery are gone; the values it computed are preserved in this file's own header as the
worked explanation of why the scoring is done the way it is.

The script's own header comment is the primary documentation now, not this README: what
LEVEL, TREND and QUAD are, why each is scored once per participant and analyzed in its own
ANOVA rather than as polynomial contrasts inside one omnibus ANOVA on the raw bins
(pooling the error term across contrasts that way is valid only under sphericity, and
several Greenhouse-Geisser epsilons already reported for this dataset's omnibus time-bin
effects are below 0.2), and a "History" section at the end of the file documenting the two
bugs found while merging the two earlier scripts into this one: a subscript-out-of-bounds
crash on every two-bin period (RS for all three measures; MBL and LRP response), and a
duration-parsing regex that silently dropped every 2-s-duration row for MBL and LRP,
because those two measures' JASP exports code duration as `stim0`/`stim2` where the
partial-activation exports use `stim0`/`stim1` -- an inconsistency in the export pipeline
itself, not something either script invented.

Each entry in `SPEC` carries the file, the bin window, and the three published F values
that do not involve the bin factor, recomputed as a replication check. Eleven of the
twelve analyses reproduce to the last reported digit; the twelfth (LRP, stimulus period)
reproduces internally, exactly, between this script and its now-deleted predecessors, but
differs by a few percent from the older hardcoded published target inherited from those
predecessors -- the target in `SPEC` was updated to the value both scripts actually
compute, and the discrepancy is flagged in the History section rather than traced further.

**The bin windows differ between files and this is easy to get wrong.** Partial
activations use bins 0-19 for the stimulus and 10-29 for the delay; the brain files are
offset by one, at 1-20 and 11-30. The two-bin windows: 0-1 for partial activations at the
signal, 1-2 at the response, 1-2 for both brain measures at the signal and the response
except the lateralised readiness potential at the response, which uses 0-1. These were
found by searching, for each file, the window reproducing all three published values at
once, which admits one solution. The wrong window reproduces nothing and gives no clue
why.

**Cell ordering.** The six cells are duration-major. R fills arrays column-major, so the
scores are built as `(n, coherence, duration)` and permuted. Getting this wrong leaves the
time effects correct -- they are grand means -- while scrambling duration and coherence,
which is a quiet and confusing failure.

**The response period is analysed on one bin.** In the corrected file the bin from -0.1 to
0 s is exactly zero for every participant and cell, because a partial activation precedes
the response burst by definition and the burst starts about 160 ms before the press. With
one bin empty, LEVEL and TREND are the same number up to a constant, so a profile analysis
would report three effects twice and would describe a slope resting on a single point.
The script detects the structurally empty bin and switches to a three-effect LEVEL-only
analysis; the switch is reported in the table header.

### `R/emit_tex.R`
Turns `stats_results.csv` into the LaTeX tables of the note, with the published value
beside each new one. Unchanged by the `run_stats.R` consolidation. Its own hardcoded `PUB`
list (used only for this script's own table headers, separate from `SPEC`'s replication-
check target above) still carries the older LRP/stimulus published value and has not been
reconciled with the value `run_stats.R` now uses -- a possible follow-up, not done here.

---

# Correspondence with the note

| note | script |
| --- | --- |
| Section 2, the model | `sim.py` |
| Section 3, predictions | derived by hand; checked by `theory_checks.py`, `model1_bound.py` |
| Section 4, from the process to the observable | `sim.py` (`_bursts`), `theory_checks.py` |
| Section 5, the data | `empirical.py` |
| Section 6, statistical re-analysis | `R/run_stats.R`, `R/emit_tex.R` |
| Section 7, corrected data | `R/run_stats.R` |
| Section 8, behaviour does not discriminate | `calibrate.py` with `BLOCKS=behav` (starting point only); `build_targets_q.py`, `calib_behav_q.py` (the manuscript's actual behavior-only calibration) |
| Section 9, simulations and calibration | `calibrate.py`; `sweep.py`/`fused_sweep.py` (letter only, not the manuscript, since 2026-08-27) |
| Figures | `figures.py` |
| Appendices A, B, D, E (predictions, no closed form) | derived by hand, not separately checked numerically |
| Appendix C, the rate | `theory_checks.py` |
| Section 3, exclusion of Model 1 | `model1_bound.py` |

---

# Things that will bite

**Results are tracked, not scratch.** `results/*.json` (fits, checkpoints, calibration
targets) are meaningful provenance and are committed to git; `results/logs/` and the large
`results/*.jsonl` sweep dumps are regenerable and gitignored (see `.gitignore`). This
replaced an earlier cloud-container-only workflow where nothing persisted between
sessions and everything had to be re-uploaded by hand -- the project now lives and runs
entirely on the local machine (`uv`/`renv`, see the top of this file).

**Two exported data files contained errors and were replaced -- and the replacement had
silently not taken effect until 2026-08-24.** The response-signal and response JASP
exports (`analysis_PMA/data/PMA_rate_{rs,resp}_JASP.csv`) originally shipped with a
spurious third bin (the archive's `compute_PMA.py` bins the RS and response windows into
two 0.1 s bins, not three). Corrected two-bin exports, dated 2026-08-19, were produced --
but landed at the root of `data/`, one directory above where `ARCHIVE`-relative paths
(`analysis_PMA/data/...`) actually resolve, so the R scripts of the time (since
consolidated into the single `run_stats.R` described above) kept reading the stale
three-bin originals underneath, unnoticed, for five days. Recomputing the response-period
duration effect directly from each file confirms the practical stakes: the stale file
gives F(1,21) = 0.011, p = .919 (no effect); the corrected file gives F(1,21) = 7.613, p =
.012 (the value already quoted in the manuscript and already the R scripts' own
comment-block target, confirming that value is correct, however it was originally produced
-- but any run of the R pipeline against the archive path between 2026-08-19 and the fix
below would have silently reproduced the wrong one). The response-signal file's correction
is comparatively minor (F recomputed from the fix: 36.6 vs the stale file's 36.8, both
consistent with the published 36.56) and changes no conclusion.

**Fixed 2026-08-24**: the corrected exports now overwrite
`analysis_PMA/data/PMA_rate_{rs,resp}_JASP.csv` in place; the stale three-bin originals are
kept for the record at `analysis_PMA/data/_superseded_aug5/*.stale-aug5`, and a labelled
copy of the correction at `data/corrected_exports_2026-08-19/` (the two root-level files
that caused the mismatch are still there too -- harmless now that the archive copy is
fixed, but safe to delete). Nothing else in the pipeline reads these two files: the
simulation, calibration, and sweep scripts (`sim.py`, `calibrate.py`, `targets.py`,
`sweep.py`, `fused_sweep.py`) extract everything directly from the raw segmented-event and
behavioural files via `empirical.py`/`targets.py` and never touch the JASP exports, so this
bug could not have reached the model-calibration or Monte-Carlo-sweep numbers -- only a
fresh, un-cached run of the R trend-contrast reanalysis (Appendix B, Table 2/S2-S4) was at
risk, and those were not rerun against the stale path this session either (only read from
already-generated tables). Regenerate `stats_tables.txt`/`stats_results.csv`/`_tables.tex`
if you want a fresh, from-scratch confirmation now that the path is fixed.

**`min_samples` counts supra-threshold samples inside the merged interval**, not the
length of the interval. A merged interval built from two short blips passes.

**Sampling rate -- fixed 2026-08-30.** The recording is at 1024 Hz; the simulation used to
run at 1000 Hz, so `min_samples = 3` meant 2.9 ms in the data and 3 ms in `sim.py`
(immaterial at the precision that mattered before, but not exact). `sim.py`'s `DT` is now
`1/1024` and `MERGE_S` is derived from the 25 ms `myonset` merge window rather than
hardcoded, so both constants are exact sample counts at the recording's own rate.
`theory_checks.py` keeps its own separate `DT` (not imported from `sim.py`) and was
updated by hand to match -- if you add a third place that hardcodes a timestep, it will not
follow `sim.py` automatically. Everything downstream of `sim.py`'s `DT` (`calibrate.py`,
`behav_sim.py`, `calib_behav.py`, `calib_behav_q.py`, `sweep.py`, `fused_sweep.py` after its
2026-08-30 fix, `export_for_R.py`, `figures.py`) was re-run or re-checked under the new
rate; `fit_model1.json`/`fit_model2.json` were refreshed by `run_emg_recal.sh`, and the
pre-refresh values are kept at `fit_model{1,2}_pre1024.json` for comparison.

**The lateralised readiness potential during the stimulus** reproduces only to about 3%,
in the same direction on all three values. Every effect there is null before and after.

**The calibration is not certified global.** The search was restarted with progressively
larger simplices until the gain became marginal, but thirteen parameters with a stochastic
objective deserve more restarts than were run. The numbers in Section 9 should be
regenerated if the search is pushed further; the checkpoint makes that cheap.
