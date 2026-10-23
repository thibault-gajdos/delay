# Analysis code

Two languages, split by job.

**Python** owns the model: the stochastic simulation, the electromyographic detector, the
parameter-space evaluation, the closed-form rate formulas, the extraction of events from
the archive, and the figures. Nothing is fitted — the models' predictions are derived
analytically and their effect-size ranges come from sampling the parameter space.

**R** owns the statistics: the analyses of variance on the electromyographic and
electroencephalographic measures, and the LEVEL, TREND and QUAD contrasts on the binned
measures. The behavioural analyses of Tables 1 and S1 are computed here too, from the
archive's own per-participant × cell exports in `analysis_behavioral_emg/` — the files
JASP was given. Accuracy and response time reproduce the published table to the last digit
(22.676 / 222.223 / 12.171 and 65.790 / 14.190 / 4.607, partial eta squared 0.519 / 0.914 /
0.367 and 0.758 / 0.403 / 0.180). Premotor and motor execution time do not quite — 69.600
against 69.75, 1.095 against 1.03 — and since the first two match exactly from the same
files by the same procedure, the difference is in the data rather than the method: the
published PMET and MET analyses appear to have run on a different version of those exports.
No verdict changes, and the manuscript's Table S1 is aligned on this output.

The two sides meet at a CSV file written by `export_for_R.py`, which carries counts and
exposures for the data.

Everything is reproducible from two inputs: the OSF archive (`86cbw`) and the manuscript.
No network access is needed once the two package managers below have fetched their
dependencies once.

## Layout

```
scripts/
  delaymodel/          the model package: sim.py, empirical.py, targets.py
  *.py                 flat, top-level entry-point scripts
  results/             the citable outputs (JSON, tracked), the sweep dumps, shards/ and
                       logs/ (all regenerable, untracked), and leak/ (the C2.3 leak reruns)
  plots/               figures written by figures.py, figures_space.py and mkfig.py
                       (untracked; the manuscript's copies are manuscript/latex/figures/)
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

`uv run python <script>` uses the locked, project-local `.venv`. `uv.lock` pins exact
versions; re-run `uv sync` after pulling a change to `pyproject.toml`. `delaymodel` is a
plain local package (no install step), resolved because every entry-point script inserts
its own directory — `scripts/`, i.e. `delaymodel`'s parent — onto `sys.path` before
importing it.

`DELAY_DATA` points at the unpacked archive. Its default is a container path, so **set it
explicitly**, e.g. in your shell profile:

```
export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
```

`SIMDIR` (R-side) points at the directory holding the CSV bridge file, default the parent
of `R/`, i.e. `scripts/`.

## R (renv)

`run_stats.R` reads the archive through `ARCHIVE`, which defaults to `../data` relative to
wherever `Rscript` is invoked from — i.e. to `scripts/data`, which does not exist. **Set it
before running**, or the first `read.csv` fails:

```
export ARCHIVE="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
```

```
cd scripts/R
Rscript -e 'renv::init()'            # once: creates renv/, renv.lock, .Rprofile
Rscript -e 'renv::install(c("afex","emmeans","dplyr","tidyr","lme4","lmerTest"))'
Rscript -e 'renv::snapshot()'        # records the installed versions in renv.lock
Rscript run_stats.R
```

On a machine that already has `renv.lock`, skip straight to `Rscript -e 'renv::restore()'`.

## Order of execution

All paths are relative to `scripts/`. Prefix any of these with `uv run python` to use the
managed virtualenv. The pipeline runs top to bottom: everything downstream of the sweep
reads its dump, and is wrapped in `./rerun.sh`, which takes about two minutes. The sweep
itself is not in it: it is long, and is run sharded across cores by `./run_sweep.sh`.

```
python delaymodel/empirical.py    -> (printout)              extraction, checked against the archive
python delaymodel/targets.py      -> (printout)              empirical summaries per cell

# 1. the parameter-space evaluation  (see fused_sweep.py below)
#    ./run_sweep.sh does this natively, one shard per core, resumable;
#    ./chunk_sweep.sh does it in bounded slices
DONE='results/shards/shard_*.jsonl' SHARD_N=$n SHARD_I=$i \
  python fused_sweep.py 1200 10000 results/shards/shard_$i.jsonl
python sweep_summary.py            -> results/sweep_screened.jsonl  the concatenated dump
                                      results/sweep_summary.json

# 1b. the other three analysis windows, for the draws the effect panel uses
#     ./chunk_windows.sh runs it in bounded slices; SHARD_N/SHARD_I split it across cores
python leak_index.py              -> results/leak/adm_idx.json     the 220 admissible draws
python augment_windows.py         -> results/window_effects.jsonl  level and trend in the
                                                                   stimulus, RS and response
                                                                   windows, at the sweep's
                                                                   own seeds

# 2. what it says, and how much of it is noise
python effect_sizes.py results/sweep_screened.jsonl 4 4000 0
                                  -> results/effect_sizes_stage1.json  distributions and
                                                                   criterion counts (the
                                                                   trailing 0 skips the
                                                                   re-simulation)
ES_SHARD_N=6 ES_SHARD_I=$i \
  python effect_sizes.py results/sweep_screened.jsonl 4 4000 40 > results/logs/es_sh$i.log
python pool_noise.py              -> results/noise_pooled.json     the variance decomposition,
                                                                   pooled over all 40 draws
                                                                   (never average the shards'
                                                                   own -- see pool_noise.py)
python ramp_check.py 3 40000 results/sweep_screened.jsonl
                                  -> results/ramp_check.jsonl      the best candidates
                                                                   re-simulated at high NTR

# 3. leak robustness (reviewer C2.3): the draws admissible at zero leak, re-run over the
#    grid gamma_x, gamma_q in {0, 0.25, 0.5, 1}. Needs the complete main sweep.
#    ./run_leak.sh does the whole grid natively; ./chunk_leak.sh does it in bounded slices
python leak_index.py              -> results/leak/adm_idx.json     the draws re-run
python leak_summary.py            -> results/leak_summary.json     verdicts at every grid
                                                                   point, and those that
                                                                   differ from zero leak

# 4. the analytic checks, which need no sampling at all
python mkfig.py                   -> plots/fig_approx_error.pdf    error of the stationary
                                                                   approximation (manuscript
                                                                   Figure 7)
python theory_checks.py           -> (printout)                    verification of the Appendix

# 5. the CSV bridge, then the figures
# order matters: figures.py reads results/representative.json, figures_space.py reads
# pma_counts.csv and results/effect_panel.json
python representative.py          -> results/representative.json   the vector Figure 2B is
                                                                   drawn at
python export_for_R.py            -> pma_counts.csv                counts and exposures, empirical
                                                                   (written at scripts/ root: R's
                                                                   SIMDIR default, "..", expects it
                                                                   there)
python figures.py                 -> plots/fig_traj_m{2,1}.pdf     manuscript Figure 1
                                     plots/fig_single.pdf          manuscript Figure 2B
python effect_panel.py            -> results/effect_panel.json     effect-size ranges vs observed
python figures_space.py           -> plots/fig_effect_sizes.pdf    those ranges, as a panel
                                     plots/fig_delay_profiles.pdf  the delay rate curve

# 6. the statistics, in R
Rscript R/run_stats.R             -> R/stats_tables.txt            all tables
                                     R/stats_results.csv           one row per effect

# 7. the manuscript's own copies, at 200 dpi. There is no script for this; the mapping is
#    fig01,02 <- fig_traj_m2,m1;  fig04 <- fig_single;  fig07 <- fig_delay_profiles;
#    fig08 <- fig_effect_sizes;  fig09 <- fig_approx_error.
#    fig05 comes from manuscript/figures/fig3.py (behaviour, drawn straight from the
#    archive); fig03, 06 and the supplemental fig10 are not produced from code here.
cd plots && pdftoppm -r 200 -png -singlefile fig_effect_sizes.pdf ../../manuscript/latex/figures/fig08
```

**The stationary start.** `sim.py` draws the maintenance noise at maintenance onset from
its own stationary law rather than setting it to zero (`WARM_START`). One line, and it is
what the delay-period analysis rests on. Started at zero, the fluctuation of the signal
about the maintained level carries a mode decaying at 2·kappa — as slow as 0.4 /s over the
sampled range — so the transient outlasts the delay and every delay-period statement is
conditional on it. Started at stationarity the transient decays at a rate of at least
lambda ≥ 10 /s, so the analysis window (which opens 1 s after stimulus offset) lies wholly
in the stationary regime, and Model 1's delay-period rate is exactly constant at every
parameter setting. The discriminating statistic is accordingly the raw-rate TREND.

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
Defines `CELLS` (the six coherence × duration conditions) and the empirical summaries per
participant and cell — accuracy, response time, premotor and motor execution time, delay
and stimulus partial-activation rates, the delay-period TREND and QUAD contrasts, the
delay log-rate slope (printed in the summary, not used downstream), premature-response
rate. These are descriptive targets used by the figures and by the
screen's observed-rate window; they are not fitted to.

TREND and QUAD come from the archive's own bin export (`analysis_PMA/data/
PMA_rate_delay_JASP.csv`) rather than from the event-file reconstruction that produces
everything else here. That export is the file the published analysis ran on, so contrasts
taken from it reproduce the manuscript's F values exactly; the reconstruction gives a mean
trend of 0.5551 against the export's 0.5604, and a 1% offset between the comparison target
and the number in the Results is not worth having.

Dispersion uses the interquartile range, not the standard deviation. The empirical
premotor execution time has a long right tail — 2.7% of trials over one second — that the
model has no mechanism for, and the untrimmed standard deviation is dominated by it.

### `sim.py`
The model of the manuscript's Theory section, and the detector.

`simulate(p, c, D, n, model, seed, gamma_x, gamma_q)` integrates the equations at 1024 Hz
(`DT = 1/1024`, the EMG/EEG native acquisition rate) and returns the gated drive of the
selected channel per trial. `model=1` sets the preparation term to zero. `gamma_x` and
`gamma_q` are the leaks of the motion-evidence accumulator and of the response-signal
accumulator, varied independently for the reviewer's sensitivity request.

`measures(sim)` returns accuracy, response time, premotor and motor execution time,
validity flags and the partial-activation onset times.

No EMG post-processing is applied to the simulation. An activation begins at every
upcrossing of the gate by the drive: no merging of nearby excursions, no minimum number
of supra-threshold samples, no refinement of onsets. The excursion that reaches the
response criterion is the response burst; its onset is where premotor time is measured
from, and it is not a partial activation. `excursions()` returns the excursions of one
trace, for Figure 2B.

### `fused_sweep.py` — the parameter-space evaluation
The evaluation of what each architecture can produce over its whole admissible parameter
space (Methods, "Evaluation over the parameter space"; Results, the "Simulation validation"
subsections). Nothing in the project is fitted, so this is where every effect-size
statement comes from.

Latin-hypercube sample over the ranges of Table 1; the production run is 1200 draws at
10000 trials per cell per model per draw (~14 s/draw on one core). At every
draw: the admissibility screen, the two RT-effect criteria, the delay-period criterion on
the raw-rate TREND (all six cells significantly positive), and the motor-execution-time
invariance criterion. The same sampled vector is evaluated in both models, with beta set
to zero in Model 1, so the two are compared at identical parameters.

Resumable by draw: a draw already recorded in `OUT`, or in any file matching the glob in
`DONE`, is skipped, so a run can be interrupted and restarted with any shard count.
Shardable: `SHARD_N`/`SHARD_I` in the environment make shard j take the draws i with
i % n == j, each writing its own file. `IDX` restricts the run to a JSON list of draws. The LHS cube is
regenerated identically from the fixed seed in every process, so `cat`-ing the shard files
reproduces the unsharded run exactly.

The module docstring carries the ranges and the provenance of every bound. Two carry most
of the weight:

  - `kappa` in [0.20, 30.0] s⁻¹, log-uniform. 1/kappa spans 5 s (slower than the whole
    delay) to 33 ms. The range is wide because nothing constrains kappa, and it favours
    neither model: the transient decays at at least lambda whatever kappa is, so no value
    of kappa gives Model 1 a delay-period ramp (Appendix, Corollary A.3).
  - `lam` in [10.0, 50.0], log-uniform. 1/lam is the time constant of the filter relating
    decision-related input to motor preparation, so the range corresponds to 100 down to
    20 ms; lam = 100 would mean 10 ms, faster than plausible neuromuscular filtering. This
    bound is a modelling commitment rather than a convenience, and is stated as one in the
    manuscript's Table 1 note: lam sets how fast the fluctuation settles, so it is what
    puts the analysis window in the stationary regime (Appendix, Proposition A.1).

No bound is derived from these data. A range drawn around a value estimated from the data
cannot test whether a prediction survives away from that value.

#### The EMG criterion in the screen

Every behavioural criterion in `admissible()` — premature rate, very fast responses,
no-responses, accuracy at 50% coherence — leaves one hole, and it is easy to miss and hard
to defend: **about half of the behaviourally admissible Model 1 draws produce exactly zero
delay-period partial activations**, and would be carried into the effect-size summaries
anyway.

The tempting fix is to reparameterise — to sample `a = (g - theta)/SDe`, the standardised
distance to the gate, instead of `SDe`, since the crossing rate goes as exp(-a²/2) and the
raw ranges put half the box above a = 3. That is the wrong instrument. Reparameterising
changes the *distribution over the parameter space*, not the acceptance set; doing it in
terms of the quantity that governs the outcome under test redefines the space in a way
chosen for its effect on the result. The screen, by contrast, acts on the acceptance set
and is neutral with respect to parameterisation.

So `SDe` is sampled directly, and `admissible()` carries one further criterion, stated at
the level where it is interpretable. A draw of parameters stands for a hypothetical
participant, so the criterion asks whether that participant's mean delay-window
partial-activation rate falls inside the range real participants span. Observed (132
participant × condition cells, none empty): 0.085–2.785 s⁻¹, mean 0.709, median 0.603;
per-participant means 0.180–1.863. The window `RATE_LO..RATE_HI = 0.04..3.0` is
deliberately wider than either, because the conclusion is an *exclusion* and erring wide is
conservative. In practice the ceiling never binds — no admissible draw of either model
exceeds 3.0 s⁻¹ — so the criterion is a floor.

Two things about the floor. **It is not a near-threshold decision.** At NTR = 4000 a draw
exactly at 0.04 s⁻¹ yields 320 events per cell and 1920 over six, a counting SE of 5.6% per
cell and 2.3% pooled; what the floor removes sits at 0 exactly, or (Model 1's ramp
candidates) at 0.001–0.005 s⁻¹, one to two orders of magnitude below it. **And it induces
no winner's curse**, because it selects on the *rate* while the reported statistic is the
*trend*: unlike a max-over-noisy-draws rule it removes the draws whose trend rests on a
handful of events, which lowers rather than inflates the noise in the reported extremes.
Any statistic that is itself a maximum over cells or draws still needs the two-stage
re-simulation `ramp_check.py` performs.

#### The chronometric criterion

A second participant-level criterion, on the LEVELS of the two components of response
time: mean premotor execution time in 181–677 ms and mean motor execution time in
46–288 ms, the observed per-participant ranges widened by a factor of two about their
centre. Without it the admissible set contains draws whose mean PMET is **−945 ms** — a
response burst beginning almost a second before the reaction signal — and whose mean MET
reaches 1.4 s. Those are not participants.

It constrains levels and never effects, which is what makes it admissible here: the
RS-locked predictions are that coherence and duration move PMET and leave MET alone, and
those remain free to fail at every admissible draw. Screening on MET *invariance* would
convert a result into an assumption, and is not done. Screening on the levels removes no
tested prediction, because the paper makes no claim about them.

It has a cost, and it should be stated. The criterion removes most of the Model 2 draws that
reached the observed MAGNITUDE of the stimulus-period duration effect. That comparison
therefore has to be made on the SIGN, which is what the derivation predicts anyway — and
there the separation is undamaged, 33% against 0%.

**Result** (`results/sweep_screened.jsonl`, 1200 draws, NTR 4000). Model 1: 160 admissible;
delay-period trend median −0.007 /s², 5th–95th percentile −0.058 to +0.004, largest value
+0.007; **0** draws with all six condition trends significantly positive and **0** anywhere
in the observed 95% confidence interval [0.375, 0.745]. Model 2: 145 admissible; trend
median +0.054, −0.017 to +0.418, largest +1.157; 78 with all six positive and 6 in the
observed interval with the rate also in range. The observed +0.5604 sits at Model 2's
96th percentile and above every draw of Model 1.

A fresh run with a different seed will not reproduce the exact admissibility counts; the
criteria and their signs are robust, the counts are not guaranteed identical run to run.
Keep the summary JSON alongside the script rather than regenerating it casually.

### `ramp_check.py` — the existence question, at high precision
An existence question is badly served by raising the trial count uniformly across a sweep:
the statistic that answers it is a MAXIMUM over draws, and a maximum over noisy estimates
is biased upward, so precision belongs on the candidates rather than spread evenly. Stage 1
is `fused_sweep.py` at moderate NTR, to find candidates; this is stage 2, re-simulating the
best of them at NTR 40000 on several independent seeds, selected two ways — the best among
those producing a realistic rate (the honest conditional question) and the best overall
regardless of rate.

The question has an exact answer — Corollary A.3: Model 1's delay-period rate is constant,
so its trend is zero everywhere — and the script checks that the simulation agrees once the
Monte-Carlo error is small enough to see the difference. It reports the verdict on exit,
recomputed from the whole file so an interrupted run still says something.

Result (`results/ramp_check.jsonl`), seven candidate draws at NTR 40000, three seeds each:

| draw | rate /s | mean trend /s² | SE |
|------|---------|-----------------|-------|
| 628  | 0.216   | +0.0007         | 0.0005 |
| 461  | 0.217   | +0.0000         | 0.0005 |
| 6    | 0.141   | −0.0002         | 0.0006 |
| 1164 | 0.120   | −0.0014         | 0.0009 |
| 15   | 0.972   | −0.0022         | 0.0002 |
| 220  | 0.216   | −0.0034         | 0.0007 |
| 965  | 1.798   | −0.0085         | 0.0028 |

The largest is +0.0007 /s² against an observed +0.5604: eight hundred times below, and
zero to within the error. The negative values are the response-burst exclusion of Appendix
Remark A.7, and they scale with the rate exactly as that account requires — the most
active draw is the most negative. This is the independent confirmation that the sweep's
+0.007 maximum is noise and not a residual ramp.

### `augment_windows.py` — the three windows outside the delay

The sweep records, per condition, a level and a linear trend for the delay window and a
level for the stimulus window. The Results report effects in the other windows too — the
coherence- and duration-dependence of the stimulus-period ramp, the fall of the rate after
the response signal and its duration-dependence. This supplies them without a full
re-sweep: only the draws admissible in either model matter to the panel, they are already
listed in `results/leak/adm_idx.json`, and each is re-simulated at the sweep's own seeds so
the new summaries belong to the same realisation as the ones in the dump. `effect_panel.py`
merges the output by draw index. `fused_sweep.window_stats` follows the empirical window
definitions exactly, including the per-bin exposure the response-signal and response
windows need, since a trial stops contributing once its response burst has begun.

One asymmetry between the two sides is worth knowing, because it decides what can be
drawn. The activation preceding a response IS the response and is excluded from the count
on both sides; in the model it occupies almost the whole 200 ms before the button press, so
the simulated pre-response rate is below 0.001 /s against an observed 0.124, and the
response-signal LEVEL is low by an order of magnitude for the same reason. Contrasts within
the response-signal window survive this (the fall of the rate after the signal is
reproduced at 94–98% of draws); levels there and everything in the response window do not,
and are not drawn.

### Leak robustness (reviewer C2.3)

The leaks of the two accumulators are varied independently: `gamma_x` on the
motion-evidence accumulator, `gamma_q` on the response-signal accumulator, over the grid
{0, 0.25, 0.5, 1} /s, every combination but (0, 0). At each point `fused_sweep.py`
re-runs, with `IDX=results/leak/adm_idx.json`, the draws admissible in either model at zero
leak, at the main sweep's trial count and seeds, into `results/leak/gx<gx>_gq<gq>/`; only
the leaks differ. Draws inadmissible at zero leak are not re-run, so a leak can remove
admissible draws but cannot add any. `leak_summary.py` classifies every grid point as
`classify.py` classifies the main sweep, and lists the verdicts (reproduced, unreproduced,
undecided) that differ from zero leak.

### `effect_sizes.py` — the citable summary, with its noise measured
Stage 1 distils the sweep dump into the counts and distributions quoted in the manuscript.
Stage 2 measures how much of the reported spread is simulation noise rather than parameter
sensitivity, by re-simulating a subset of admissible draws on several seeds and decomposing
the variance of the per-draw trend into between-draw and within-draw (Monte-Carlo) parts.

The statistic throughout is TREND — the linear orthogonal polynomial contrast of the raw
rate, in /s². A log-rate slope is a ratio, so a draw can earn a steep one by ramping from
almost nothing; the trend is the product of the level and the relative rise, so it demands
both, which is what the observation shows. Being linear in the rate it also has one
unambiguous observed value, +0.5604 /s², whether taken as the mean of the per-cell,
per-participant contrasts or as the contrast of the pooled profile.

The acceptance band is the between-participant 95% confidence interval on the observed
trend, [0.375, 0.745] /s², and not the spread of the six condition means ([0.469, 0.634],
still reported alongside): a spread across conditions is not an uncertainty on the quantity
a model is asked to reproduce. Both counts are in the output, as are two prior-free
summaries — `pct_of_observed` (Model 2: 96) and `ratio_observed_to_max` (Model 1: 75) —
which are what the Results lead with, because a count of draws in a band reports the
Latin-hypercube prior as much as the architecture.

Shardable through `ES_SHARD_N`/`ES_SHARD_I`. Pool the per-draw values from the shard logs
with `pool_noise.py` rather than averaging the shards' own decompositions: a between-draw
variance estimated from seven draws is worthless, and one shard duly reported a
Monte-Carlo share above 100%.

Result: Model 1 total SD 0.099 s⁻¹, of which 0.055 Monte-Carlo and 0.083 between draws
(31% of the variance); Model 2 total 1.18, of which 0.057 Monte-Carlo and 1.17 between
draws (0.2%). Model 1's trend is not merely centred near zero but pinned there — the whole
span of its admissible space moves it about a tenth as much as Model 2's space moves
Model 2's.

### `effect_panel.py` — what each architecture can produce, on every reported effect
This is what the parameter-space evaluation is FOR, and the manuscript's Figure 6. It
covers every effect on behavior or partial-activation rate that the Results report — not
only the ones that discriminate — grouped by analysis window in the order the Results take
them (`GROUPS` in the script drives both the JSON and the figure). The reason is
navigational: a reader checking a sentence of the Results should find that sentence's
effect in the figure rather than discover it was left out for being undiscriminating.

Three classes are absent, and the figure caption says so: the electroencephalographic
measures, for which the models predict a direction but no magnitude; a few interactions
null in the data and in both models; and the response window together with the LEVEL in the
response-signal window, for the reason given under `augment_windows.py`.

The fraction of admissible settings carrying the OBSERVED SIGN is the compact summary of the
whole comparison:

| effect | Model 1 | Model 2 |
| --- | --- | --- |
| the effects both architectures share | 86–98% | 88–94% |
| **duration effect on the stimulus-period rate** | **0%** | **33%** |
| **duration effect on the stimulus-period ramp** | **0%** | **37%** |
| **the delay-period rate trend** | **1%** | **80%** |
| coherence × duration on the delay rate | 3% | 5% |

Three readings worth keeping. **The stimulus window excludes Model 1 twice**, on the level
and on the ramp, by sign: its whole admissible range is negative where the observation is
positive, and both are independent of everything in the delay period. Note that the test is
ONE-SIDED — temporal preparation opposes the level effect but overturns it at only a third
of Model 2's draws, because the ramp's size and the distance to the gate are independent
parameters — so it refutes Model 1 and does not confirm Model 2. **The coherence × duration
interaction on the delay rate is missed by both models**, across their whole admissible
spaces: the paper's acknowledged disagreement, established as a property of the
architectures rather than of a parameter choice. And the **PMET/MET dissociation is
confirmed in the data but is a weak test**: both architectures can also produce
motor-component effects far larger than anything observed, so finding a small one does not
discriminate.

### `representative.py` — the point the illustrative figures are drawn at
Nothing is fitted, so the figures that show simulated dynamics need a rule for choosing one
parameter vector. It is the admissible draw closest to the coordinate-wise median of the
sampled parameters of the draws admissible in BOTH models, each coordinate standardised by
its interquartile range and taken in logarithm where it was sampled in logarithm. The rule
touches only the parameters — never the rate, the trend, or anything the figures are used
to discuss.

Two traps, both recorded in the module docstring. The coordinate-wise median of the
*ranges* is not usable: it sits at a standardised distance to the gate near 4, where the
model produces no muscle activity at all. And beta must be set separately, to the median of
the values admissible under Model 2, because requiring admissibility in both models biases
beta downward (0.0059 against 0.0088) — a large ramp is what produces premature responses,
and the screen removes those — so taking beta from the draw would illustrate Model 2 at a
ramp smaller than is typical of Model 2, for a reason that has nothing to do with Model 2.

### `figures_space.py`
`fig_effect_sizes.pdf` is the panel described under `effect_panel.py`, grouped by analysis
window. `fig_delay_profiles.pdf` shows the delay-period rate curve itself, each draw divided
by its own mean over the window so that draws whose rates differ by two orders of magnitude
are comparable on shape.

### `export_for_R.py`
Writes the bridge file. `pma_counts.csv` has one row per participant, cell, period and
bin, with the raw count and the number of trials contributing to that bin.

Exposure is not always the number of valid trials. A partial activation precedes the
response burst by definition, so once that burst has begun a trial can no longer
contribute and later bins are undefined rather than empty. This matters only in the
response-signal and response windows.

### `figures.py`
The figures that show simulated dynamics: the manuscript's Figure 1, drawn for each model
over thirty admissible draws — the same thirty in both, differing only in beta — and
Figure 2B, one simulated trial with the detected intervals, at the representative vector.
Run `representative.py` first (Figure 1 does not need it; Figure 2B does).

### `mkfig.py`
The manuscript's Figure 7: how far the stationary rate is from the exact rate, over the
admissible draws and the whole delay window. Analytic, no simulation. `s2` is each draw's
own SD_eps² — it does not cancel from the ratio, because it sets sigma_H and hence the
standardized distance to the gate, and the rate depends on that distance through a
Gaussian, so the error is amplified by a factor of order a² (Appendix, Corollary A.2).

### `theory_checks.py`
Verification of the Appendix against direct integration at a 1024 Hz step (its own local
`DT`, independent of `sim.py`'s — keep the two in sync by hand): the Kac–Rice expression
against a simulated count of upcrossings, the occupation identity against the simulated
time above the level, the ratio of merged clusters to raw upcrossings, and the residual
condition dependence of motor execution time from the filter transient.

The Euler bias on crossing counts is upward, about 2% at a 1 ms step and shrinking with
the step. It affects the simulated rate, not the direction of any effect.

# R

### `R/run_stats.R`
Every test: the LEVEL / TREND / QUAD decomposition, written with `afex` (repeated-measures
ANOVA, Greenhouse–Geisser) and `emmeans` (trend contrasts). Writes `stats_tables.txt` and
`stats_results.csv`. It analyses the empirical counts only: a parameter vector is one
hypothetical participant, so the simulated design has no between-participant factor, and
what the models predict is reported as effect-size ranges over the parameter space.

The script's own header is the primary documentation, not this README: what LEVEL, TREND
and QUAD are, and why each is scored once per participant and analyzed in its own ANOVA
rather than as polynomial contrasts inside one omnibus ANOVA on the raw bins — pooling the
error term across contrasts that way is valid only under sphericity, and several
Greenhouse–Geisser epsilons for this dataset's omnibus time-bin effects are below 0.2.

Each entry in `SPEC` carries the file, the bin window, and the three published F values
that do not involve the bin factor, recomputed as a replication check. Eleven of the twelve
analyses reproduce to the last reported digit; the twelfth (LRP, stimulus period) is
discussed at the end of the script.

**The bin windows differ between files and this is easy to get wrong.** Partial
activations use bins 0–19 for the stimulus and 10–29 for the delay; the brain files are
offset by one, at 1–20 and 11–30. The two-bin windows: 0–1 for partial activations at the
signal, 1–2 at the response, 1–2 for both brain measures at the signal and the response
except the lateralised readiness potential at the response, which uses 0–1. These were
found by searching, for each file, the window reproducing all three published values at
once, which admits one solution. The wrong window reproduces nothing and gives no clue
why.

**Cell ordering.** The six cells are duration-major. R fills arrays column-major, so the
scores are built as `(n, coherence, duration)` and permuted. Getting this wrong leaves the
time effects correct — they are grand means — while scrambling duration and coherence,
which is a quiet and confusing failure.

**The response period is analysed on one bin.** The bin from −0.1 to 0 s is exactly zero
for every participant and cell, because a partial activation precedes the response burst by
definition and the burst starts about 160 ms before the press. With one bin empty, LEVEL
and TREND are the same number up to a constant, so a profile analysis would report three
effects twice and would describe a slope resting on a single point. The script detects the
structurally empty bin and switches to a three-effect LEVEL-only analysis; the switch is
reported in the table header.

### `R/emit_tex.R`
Turns `stats_results.csv` into LaTeX tables, with the published value beside each new one.
Its own hardcoded `PUB` list (used only for this script's own table headers, separate from
`SPEC`'s replication-check target) carries the older LRP/stimulus published value and has
not been reconciled with the value `run_stats.R` computes.

---

# Correspondence with the manuscript

| manuscript | script |
| --- | --- |
| Theory (Materials and Methods) | `sim.py` |
| Theoretical predictions, each window | derived by hand; checked by `theory_checks.py` |
| Appendix, Proposition A.1 to Remark A.7 | derived by hand; `theory_checks.py`, `mkfig.py` |
| Exclusion of Model 1 | Corollary A.3 (analytic); `fused_sweep.py` + `ramp_check.py` (sampled) |
| From the process to the observable | `sim.py` (`measures`), `theory_checks.py` |
| The data | `empirical.py`, `targets.py` |
| Results, all statistics | `R/run_stats.R`, `R/emit_tex.R` |
| Evaluation over the parameter space | `fused_sweep.py`, `augment_windows.py`, `effect_sizes.py`, `effect_panel.py`, `ramp_check.py` |
| Leak robustness (letter, C2.3) | `leak_index.py`, `fused_sweep.py`, `leak_summary.py` |
| Figures | `figures.py`, `figures_space.py`, `mkfig.py` |

---

# Things that will bite

**Results are tracked, not scratch.** `results/*.json` (the effect-size summary, the effect
panel, the noise decomposition, the leak summary, the representative vector) are the
citable numbers and are committed to git; `results/logs/`, `results/shards/`, `plots/` and
the large `results/*.jsonl` dumps are regenerable and gitignored (see `.gitignore`), except
`ramp_check.jsonl`, which is small and is the record of the high-precision re-simulation.

**Two of the archive's exported data files were wrong, and the corrected copies are what
the archive now holds.** The response-signal and response JASP exports
(`analysis_PMA/data/PMA_rate_{rs,resp}_JASP.csv`) originally shipped with a spurious third
bin — the archive's `compute_PMA.py` bins those windows into two 0.1 s bins, not three —
and the response-locked file was referenced to the onset of the response burst rather than
to the button press. The corrected two-bin exports of 2026-08-19 are in place, and the
originals are kept beside them at `analysis_PMA/data/_superseded_aug5/*.stale-aug5`. It
matters: the stale response file gives F(1,21) = 0.011, p = .919 for the duration effect,
the corrected one F(1,21) = 7.613, p = .012, which is what the manuscript reports. Nothing
else in the pipeline reads these two files — the Python side extracts everything from the
raw segmented-event and behavioural files — so only the R reanalysis was ever at risk.
**Check that the deposited OSF archive carries the corrected copies**; the response letter
tells the reviewer it does.

**One timestep, two places.** `sim.py`'s `DT` is `1/1024`, the recording's own rate.
`theory_checks.py` keeps its own separate `DT` and must
be kept in step by hand; a third place that hardcodes a timestep will not follow `sim.py`
automatically either.

**The lateralised readiness potential during the stimulus** reproduces only to about 3%,
in the same direction on all three values. Every effect there is null before and after.

**The curvature is an exact prediction that the simulation cannot measure cleanly, and the
manuscript says so.** The two models make opposite and exact claims about the quadratic
contrast of the delay-period rate: zero under Model 1, since a constant rate has zero
curvature, and positive under Model 2 wherever the gate is more than one fluctuation
standard deviation away. The data give +0.2099 /s³, F(1,21) = 18.891. The simulated
curvatures do not reach it: 0 of 160 Model 1 draws and 19 of 145 Model 2 draws have all six
condition curvatures significantly positive, and both medians sit near zero.

The reason is measurement, not model, and it is worth recording because it is easy to
misread as a failed prediction. `measures()` and `empirical.py` both score the burst
preceding the response as the response rather than as a partial activation, so a trial
whose final approach begins before the response signal loses an event from the END of the
window — and the quadratic contrast weights the ends. Restricted to the first eighteen
bins Model 2's median curvature turns positive and Model 1's halves. The Results therefore
report the curvature as the analytic prediction it is and the simulated distribution as a
lower bound on it.

**Do not pool counts across participants before taking a curvature.** Pooling flattens it
on any log scale (log-sum-exp of concave functions is not concave — measured, −0.095 per
participant against −0.034 pooled). The raw-rate contrasts used here do not have this
problem, being linear in the rate, but the trap is one step away.
