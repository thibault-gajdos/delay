# Analysis code for the theoretical note

Two languages, split by job.

**Python** owns the model: the stochastic simulation, the electromyographic detector, the
parameter-space evaluation, the closed-form frontier, the extraction of events from the
archive, and the figures. Nothing is fitted -- the models' predictions are derived
analytically and their effect-size ranges come from sampling the parameter space, not from
estimating a parameter vector.

**R** owns the statistics that this project computes: the analyses of variance on the
electromyographic and electroencephalographic measures, the trend decompositions, and the
Poisson analyses of the partial-activation counts.

The behavioural analyses of variance of Tables 1 and S1 -- accuracy, response time, premotor
and motor execution time -- were until 2026-09-01 the one set of analyses in the paper that
no script here reproduced: they came from the original JASP runs and lived nowhere else.
They are now computed by `run_stats.R` from the archive's own per-participant x cell
exports in `analysis_behavioral_emg/`, which are the files JASP was given. Accuracy and response time reproduce the published table to the last digit (22.676 /
222.223 / 12.171 and 65.790 / 14.190 / 4.607, with partial eta squared 0.519 / 0.914 /
0.367 and 0.758 / 0.403 / 0.180). Premotor and motor execution time do not quite -- 69.600
against 69.75 and 1.095 against 1.03 -- and since the first two match exactly from the same
files by the same procedure, the difference is in the data rather than in the method: the
published PMET and MET analyses presumably ran on an earlier version of those exports. No
verdict changes. **The manuscript's Table S1 and the corresponding sentences of the Results
have been aligned on this output**, which also fills in the partial eta squared the table
left blank for MET (0.050, 0.029, 0.104) and replaces the uncorrected interaction degrees
of freedom (2, 42) with the Greenhouse-Geisser ones afex reports.

Two bugs found on the way, both worth knowing about:

* afex's default effect size is GENERALISED eta squared, and this file was collecting it
  into a column named `eta2p`. That is a mislabelling, not a rounding difference: for the
  accuracy duration effect the two are 0.027 and 0.519, and it is the second that the
  manuscript reports. Fixed with `es_aov = "pes"`; every row of `stats_results.csv` was
  affected, not only the new ones.
* `analyse_period` crashed on a period with NO events in any bin: `which(!empty)[1]` returns
  NA, the filter kept no rows, and aov_ez failed with an opaque reshape error. The case is
  not pathological -- the SIMULATED response period contains no partial activations by
  construction, since `measures()` scores the burst preceding the response as the response.
  It had never fired because `--with-sim` had never run to completion. The two sides meet at CSV files
written by `export_for_R.py`, which carry counts and exposures for the data and for each
model simulated at the representative parameter vector, so that the identical statistical
analysis runs on both sides.

Everything is reproducible from two inputs: the OSF archive (`86cbw`) and the manuscript.
No network access is needed once the two package managers below have fetched their
dependencies once.

## Layout

```
scripts/
  delaymodel/          the model package: sim.py, empirical.py, targets.py
  *.py                 flat, top-level entry-point scripts
  results/             the citable outputs (JSON/JSONL, tracked) and logs/ (regenerable)
  plots/               figures written by figures.py and figures_space.py
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

`run_stats.R` reads the archive through `ARCHIVE`, which defaults to `../data` relative to
wherever `Rscript` is invoked from -- i.e. to `scripts/data`, which does not exist. **Set it
before running**, or the first `read.csv` fails with "impossible d'ouvrir la connexion":

```
export ARCHIVE="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
```

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
the managed virtualenv. The pipeline runs top to bottom: everything downstream of the sweep
reads its dump.

Everything downstream of the sweep and the frontier is wrapped in `./rerun.sh`, which runs
the steps below in order and takes about two minutes. The sweep and the frontier are not in
it: both are long and both are run sharded across cores by hand.

```
python delaymodel/empirical.py    -> (printout)              extraction, checked against the archive
python delaymodel/targets.py      -> (printout)              empirical summaries per cell

# 1. the parameter-space evaluation  (~12 min over six shards; see fused_sweep.py below)
SHARD_N=6 SHARD_I=$i \
  python fused_sweep.py 600 4000 results/sweep_screened_sh$i.jsonl
cat results/sweep_screened_sh*.jsonl > results/sweep_screened.jsonl

# 2. what it says, and how much of it is noise
python effect_sizes.py results/sweep_screened.jsonl 4 4000 0
                                  -> results/effect_sizes_stage1.json  distributions and
                                                                   criterion counts (stage 1
                                                                   only; the trailing 0 skips
                                                                   the re-simulation)
ES_SHARD_N=6 ES_SHARD_I=$i \
  python effect_sizes.py results/sweep_screened.jsonl 4 4000 40 > results/logs/es_sh$i.log
python pool_noise.py              -> results/noise_pooled.json     the variance decomposition,
                                                                   pooled over all 40 draws
                                                                   (never average the shards'
                                                                   own -- see pool_noise.py)
python ramp_check.py 3 40000 results/sweep_screened.jsonl
                                  -> results/ramp_check.jsonl      the best candidates
                                                                   re-simulated at high NTR

# 3. the analytic frontier, which needs no sampling at all
python frontier.py                -> results/frontier.json         Model 1's ramp-rate frontier
python theory_checks.py           -> (printout)                    verification of Appendix C

# 4. figures
python representative.py          -> results/representative.json   the point the illustrative
                                                                   figures are drawn at
python figures.py                 -> plots/fig_traj_m{1,2}.pdf     simulated dynamics
                                     plots/fig_single.pdf
                                     plots/fig_rates.pdf
                                     plots/fig_shape.pdf
python effect_panel.py            -> results/effect_panel.json     effect-size ranges vs observed
python figures_space.py           -> plots/fig_effect_sizes.pdf    those ranges, as a panel
                                     plots/fig_frontier.pdf        the frontier with the draws
                                     plots/fig_delay_profile.pdf   the delay rate curve

# 5. the statistics, in R
python export_for_R.py            -> pma_counts.csv                counts and exposures, empirical
                                     sim_counts_mM.csv             the same from each model at the
                                                                   representative vector (written at
                                                                   scripts/ root: R's SIMDIR default,
                                                                   "..", expects them there)
Rscript R/run_stats.R --with-sim  -> R/stats_tables.txt            all tables, data and simulation
                                     R/stats_results.csv           one row per effect
                                                                   (includes the behavioural
                                                                   ANOVAs of Tables 1 and S1
                                                                   since 2026-09-01)
Rscript R/emit_tex.R              -> _tables.tex                   the same, for the note
```

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
Defines `CELLS` (the six coherence x duration conditions) and the empirical summaries per
participant and cell -- accuracy, response time, premotor and motor execution time, delay
and stimulus partial-activation rates, delay log-rate slope, premature-response rate. These
are descriptive targets used by the figures and by the screen's observed-rate window; they
are not fitted to.

Dispersion uses the interquartile range, not the standard deviation. The empirical
premotor execution time has a long right tail -- 2.7% of trials over one second -- that the
model has no mechanism for, and the untrimmed standard deviation is dominated by it.

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

### `fused_sweep.py` -- the parameter-space evaluation
The manuscript's evaluation of what each architecture can produce over its whole
admissible parameter space (Methods, "Evaluation over the parameter space"; Results,
Delay-locked). Nothing in the project is fitted, so this is where every effect-size
statement comes from.

Latin-hypercube sample over the ranges of Table A1; production runs use 600 draws at 4000
trials per cell per model per draw (~6.5 s/draw, ~12 min over six shards). At every draw:
the admissibility screen, the two RT-effect criteria (coherence, duration), the
delay-period log-rate slope criterion (all six cells significantly positive) and the
motor-execution-time invariance criterion. The same sampled vector is evaluated in both
models, with beta set to zero in Model 1, so the two are compared at identical parameters.

Resumable by line count of `OUT`. Shardable: `SHARD_N`/`SHARD_I` in the environment make
shard j take the draws i with i % n == j, each writing its own file. The LHS cube is
regenerated identically from the fixed seed in every process, so `cat`-ing the shard files
reproduces the unsharded run exactly -- verified by re-running the sweep sharded and
recovering identical per-draw statistics.

```
for i in 0 1 2 3 4 5; do
  SHARD_N=6 SHARD_I=$i uv run python3 fused_sweep.py 600 4000 \
    results/sweep_screened_sh$i.jsonl > results/logs/sweep_sh$i.log 2>&1 &
done
cat results/sweep_screened_sh*.jsonl > results/sweep_screened.jsonl
```

The optional fourth argument, `GAMMA`, threads a common leak into every `simulate()` call
(`sim.py`'s `gamma`) for the leak-robustness check of the response letter (C2.3). It has
not been rerun at the current ranges.

**Why the ranges are what they are.** Every bound is either a fixed architectural constant
(`sigx`, `Lmotion`, `LRS`, `delay`, and `Lmotor` at Dendauw et al. 2024's 0.093 s) or is
derived from the design, from a published estimate, or from the structure of the model. No
bound is derived from these data -- a range drawn around a value estimated from the data
cannot test whether a prediction survives away from that value. Two bounds carry most of
the weight and deserve to be visible:

  - `kappa` in [0.20, 30.0] s^-1, log-uniform. 1/kappa spans 5 s (slower than the whole
    delay) to 33 ms. The range is deliberately generous, and generous *towards Model 1*:
    small kappa is Model 1's only route to a delay-period rate still rising at the end of
    the delay, so narrowing it would delete the region most favourable to the model the
    evaluation exists to exclude. A floor at kappa = 2.0 was briefly imposed on 2026-08-30
    and reverted the next day; the justification offered for it -- that the stationary
    crossing-rate approximation loses validity below it -- was not a reason to exclude the
    region, because nothing in the paper uses that approximation there (`kernel_moments.py`
    computes from the kernel). For a range whose purpose is exclusion, erring wide is the
    conservative direction.
  - `lam` in [10.0, 50.0], log-uniform. 1/lam is the time constant of the filter relating
    decision-related input to motor preparation, so the range corresponds to 100 down to
    20 ms; lam = 100 would mean 10 ms, faster than plausible neuromuscular filtering. This
    bound is load-bearing in the other direction: Model 1's ramp-rate frontier
    (`frontier.py`) is most favourable at the largest lam, so capping lam is what makes
    the frontier a bound rather than an artefact of where the range stops. It is stated in
    the manuscript's Table A1 note as a modelling commitment for exactly that reason.

#### The EMG criterion in the screen

Every criterion in `admissible()` was originally behavioural: premature rate, very fast
responses, no-responses, accuracy at 50% coherence. None mentioned muscle activity, and
the consequence was easy to miss and hard to defend: **nearly half of the behaviourally
admissible Model 1 draws produced exactly zero delay-period partial activations**, and
were carried into the effect-size summaries anyway.

The near-miss fix was to reparameterise -- to sample `a = (g - theta)/SDe`, the
standardised distance to the gate, instead of `SDe`, since the crossing rate goes as
exp(-a^2/2) and the raw ranges put half the box above a = 3. That was implemented and then
rejected, correctly. Reparameterising changes the *distribution over the parameter space*,
not the acceptance set; doing it in terms of the quantity that governs the outcome under
test redefines the space in a way chosen for its effect on the result, which is a weak form
of the circularity the range revision exists to remove. The screen, by contrast, acts on
the acceptance set and is neutral with respect to parameterisation.

So `SDe` is sampled directly, and `admissible()` carries one further criterion, stated at
the level where it is interpretable. A draw of parameters stands for a hypothetical
participant, so the criterion asks whether that participant's mean delay-window
partial-activation rate falls inside the range real participants span. Observed (132
participant x condition cells, none empty): 0.085-2.785 s^-1, mean 0.709, median 0.603;
per-participant means 0.180-1.863. The window `RATE_LO..RATE_HI = 0.04..3.0` is
deliberately wider than either, because the conclusion is an *exclusion* and erring wide is
conservative. In practice the ceiling never binds -- no admissible draw of either model
exceeds 3.0 s^-1 -- so the criterion is a floor.

Two things about the floor. **It is not a near-threshold decision.** At NTR = 4000 a draw
exactly at 0.04 s^-1 yields 320 events per cell and 1920 over six, a Poisson SE of 5.6% per
cell and 2.3% pooled; what the floor removes sits at 0 exactly, or (Model 1's ramp
candidates) at 0.001-0.005 s^-1, one to two orders of magnitude below it. **And it induces
no winner's curse**, because it selects on the *rate* while the reported statistic is the
*slope*: unlike a max-over-noisy-draws rule it removes the draws whose slope rests on a
handful of events, which lowers rather than inflates the noise in the reported extremes.
Any statistic that is itself a maximum over cells or draws still needs the two-stage
re-simulation `ramp_check.py` performs.

Cost: about 40% fewer usable draws per 250 than the rejected reparameterisation, bought
back by raising `N_DRAWS` to 600.

#### The chronometric criterion (2026-09-01)

The screen gained a second participant-level criterion, on the LEVELS of the two components
of response time: mean premotor execution time in 181-677 ms and mean motor execution time
in 46-288 ms, the observed per-participant ranges widened by a factor of two about their
centre. Without it the admissible set contained draws whose mean PMET was **-945 ms** -- a
response burst beginning almost a second before the reaction signal -- and whose mean MET
reached 1.4 s. Those are not participants.

It constrains levels and never effects, which is what makes it admissible here: the
RS-locked predictions are that coherence and duration move PMET and leave MET alone, and
those remain free to fail at every admissible draw. Screening on MET INVARIANCE was
considered and rejected for exactly that reason -- it would have converted a result into an
assumption. Screening on the levels removes no tested prediction, because the paper makes no
claim about them.

It has a cost, and it should be stated. The criterion removes most of the Model 2 draws that
reached the observed MAGNITUDE of the stimulus-period duration effect: 23 of 169 before, 2
of 148 after. The B.2 comparison therefore has to be made on the SIGN, which is what
Proposition B.2 predicts anyway -- and there the separation is undamaged, 44% against 1%.

**Result** (`results/sweep_screened.jsonl`, 1200 draws, NTR 4000). Model 1: 164 admissible;
delay slope median -0.02, 5th-95th percentile -0.20 to +0.09; the observed sign in 31% of
settings; **0** draws with all six cell slopes significantly positive and **0** anywhere in
the observed slope band [0.659, 1.029]. Model 2: 148 admissible; slope median +0.27, -0.06
to +1.98; the observed sign in 90%; 76 with all six positive and 16 in the observed band.

The fraction of admissible settings carrying the OBSERVED SIGN is the compact summary of the
whole comparison, and `effect_panel.py` now reports it:

| effect | Model 1 | Model 2 |
| --- | --- | --- |
| the seven effects both architectures share (A.2, A.3, A.4, A.5, B.1, D.2 x2) | 96-99% | 96-99% |
| **B.2, duration effect on the stimulus-period rate** | **1%** | **45%** |
| **C.2, delay-period log-rate slope** | **31%** | **90%** |
| C.3, coherence x duration on the delay rate | 5% | 9% |

Two effects separate the architectures, one is missed by both, seven are indistinguishable.

**Proposition B.2 was half wrong and has been corrected.** It claimed that Model 2 predicts a
higher stimulus-locked rate after the SHORT stimulus and Model 1 the opposite. Model 1's half
holds at 99% of its admissible settings. Model 2's does not: temporal preparation opposes the
level effect but overturns it at only 45%, because the ramp's size and the distance to the
gate are independent parameters. The test is therefore ONE-SIDED -- it refutes Model 1 and
does not confirm Model 2 -- and the proposition, the Results paragraph and Remark B.3 now say
so. This was caught by reading the effect panel, which showed both models' medians on the
same side of zero where the text said they should be on opposite sides.

A fresh run with a different seed will not reproduce the exact admissibility counts; the
criteria and their signs are robust, the counts are not guaranteed identical run to run.
Keep the summary JSON alongside the script rather than regenerating it casually.

### `kernel_moments.py`
The exact non-stationary second moments of the maintained fluctuation V, and the Kac-Rice
upcrossing rate that follows from them, keeping the V/Vdot covariance that non-stationarity
induces. Computed from the kernel rather than from the stationary limits because the
process is not stationary over most of the sampled range of kappa: the relative departure
one second into maintenance is of order exp(-2 kappa), below 1e-9 only above kappa ~ 11,
whereas Table A1 samples down to 0.2. Small kappa is also the region most favourable to
Model 1, so it is exactly the region a bound must not approximate away.

(Renamed from model1_bound.py on 2026-08-31 -- that file also carried a ceiling argument
and two constants taken from the discarded fit; the ceiling is superseded by the frontier
below. The old name survives only in the git history.)

### `frontier.py` -- the ramp-rate frontier of both models
Sharded over the lambda grid: `FR_MODEL`, `FR_SHARD_N`, `FR_SHARD_I` write
`results/frontier_m{model}_sh{i}.npz`, and `python3 frontier.py merge` combines them into
`results/frontier.json`. Sharding by lambda is exact -- the scatter is a product grid, no
shard needs another's points, and the envelope of a union is the envelope of the
concatenation.
The analytic core of the delay-period argument, and the strongest thing in the project,
because it involves no sampling at all.

With alpha = 0 the delay-period log-rate slope is

    d(ln nu)/dtau = 0.5*(u'/u - s'/s) + (a^2/2)*(s'/s),    s = Var V, u = Var Vdot,

whose first term is NEGATIVE throughout the analysis window -- u is governed by the fast
filter and has saturated, s carries the slow mode and has not -- so the ramp lives entirely
in the second term, which needs a large a, while the rate carries exp(-a^2/2), which needs
a small one. **Ramp and rate are governed by the same quantity in opposite directions**,
and the ceiling omega/2pi is attained at a = 0, where the slope is strictly negative: the
configuration maximising Model 1's rate is one in which that rate declines.

The script computes the frontier as the upper envelope of achievable (rate, slope) pairs
over a grid in (kappa, lambda, a) -- faster than a constrained optimisation at each rate,
and it returns the (kappa, lambda, a) attaining each point. Queries are answered from the
raw scatter rather than from the binned envelope: binning can only lose the maximum, and a
frontier used to exclude a model must not be under-reported.

Model 2 is computed the same way with the barrier moving: temporal preparation enters the
Kac-Rice formula only through the relative derivative Vdot + alpha, so `rate(..., alpha=)`
handles it. Its ramp height is capped at 3.0 sigma_V -- the 75th percentile of beta/sigma_V
over the Model 2 draws the screen actually accepts (median 1.0, max 31.2) -- because the
analytic calculation knows nothing about admissibility, and without a cap Model 2's
frontier grows without bound simply by making the ramp larger.

Result (`results/frontier.json`), against an observed slope of +0.892 /s at a rate of
0.709 /s, with the slope estimated as the data's slope is estimated (an offset-Poisson
regression over the same twenty bins, not a secant -- see below):

| rate at least | Model 1 max slope | Model 2 max slope |
| --- | --- | --- |
| 0.085 /s | +0.81 | +5.16 |
| 0.200 | +0.40 | +4.73 |
| 0.400 | +0.05 | +4.35 |
| **0.709** | **+0.000** | **+4.05** |
| 1.000 | +0.000 | +3.83 |

and the largest rate at which the observed slope is attainable at all: **0.069 /s for
Model 1, 6.44 /s for Model 2**. Model 1's figure is the sharper statement, because 0.069 /s
is BELOW the lowest rate recorded in any participant x condition cell (0.085): there is no
admissible (kappa, lambda, a) at which Model 1 produces both the observed ramp and partial
activations as often as the least active participant in the least active condition. Model 2
at the observed rate exceeds the observed slope more than fourfold.

**The slope statistic matters, and using the wrong one was generous to Model 1.** Until
2026-09-01 the frontier reported the secant between the ends of the window while the data
were analysed with an offset-Poisson regression over twenty 0.1 s bins. For Model 1's
strongly concave log-rate the regression weights the late, event-rich bins -- where the
variance transient has already flattened -- and returns 45-75% of the secant. Correcting it
tightened every entry above (at 0.085, +1.17 -> +0.81) and moved the largest rate at the
observed slope from 0.140 to 0.069 /s, i.e. from "a fifth of the observed mean" to "outside
the observed range altogether".

The two frontiers have different STATUS, and the difference should be stated wherever they
are quoted together. For Model 1 the frontier is a bound: with alpha = 0 the rate depends
on nothing but kappa, lambda and the distance to the gate, all three swept, and nothing
outside the calculation can raise it. For Model 2 it is an upper envelope under a cap read
off the screened sample -- what the architecture can reach, not what it typically does.

### `ramp_check.py` -- the existence question, at high precision
An existence question is badly served by raising the trial count uniformly across a sweep:
the statistic that answers it is a MAXIMUM over draws, and a maximum over noisy estimates
is biased upward, so precision belongs on the candidates rather than spread evenly. Stage 1
is `fused_sweep.py` at moderate NTR, to find candidates; this is stage 2, re-simulating the
best of them at NTR 40000 on several independent seeds, selected two ways -- the best among
those producing a realistic rate (the honest conditional question) and the best overall
regardless of rate.

Result (`results/ramp_check.jsonl`): Model 1's draws at realistic rates (0.19-1.10 s^-1)
give slopes between -0.007 and 0.000 with MC standard errors of 0.004 or less -- not a
small ramp, no ramp. Its draws that do ramp (+0.10, +0.13, +0.24) do so at 0.040-0.051
s^-1. Two draws whose first-pass slopes looked positive (+0.18, +0.22) fell to +0.01 and
-0.03 on re-simulation, which is what the upward bias of a maximum looks like when it is
checked. The sampled space and the analytic frontier agree by independent routes.

### `effect_sizes.py` -- the citable summary, with its noise measured
Stage 1 distils the sweep dump into the counts and distributions quoted in the manuscript.
Stage 2 measures how much of the reported spread is simulation noise rather than parameter
sensitivity, by re-simulating a subset of admissible draws on several seeds and decomposing
the variance of the per-draw slope into between-draw and within-draw (Monte-Carlo) parts.
Shardable through `ES_SHARD_N`/`ES_SHARD_I`. Pool the per-draw values from the shard logs
with `pool_noise.py` rather than averaging the shards' own decompositions: a between-draw
variance estimated from seven draws is worthless, and one shard duly reported a
Monte-Carlo share above 100%.

Result: Model 1 total SD 0.099 s^-1, of which 0.055 Monte-Carlo and 0.083 between draws
(31% of the variance); Model 2 total 1.18, of which 0.057 Monte-Carlo and 1.17 between
draws (0.2%). Model 1's slope is not merely centred near zero but pinned there -- the whole
span of its admissible space moves it about a tenth as much as Model 2's space moves
Model 2's.

### `effect_panel.py` -- what each architecture can produce, on the effects that discriminate
This is what the parameter-space evaluation is FOR, and the paper's Figure S6. The effects
are chosen to be the contrasts that carry an argument -- the interactions, the one
prediction outside the delay on which the models disagree, and the PMET/MET dissociation --
rather than main effects both architectures produce with the same sign everywhere. An
earlier version reported accuracy at 50% coherence; it is gone, because the admissibility
screen requires accuracy above 0.80 in those cells and so pins both models near ceiling
there, which tests nothing.

Result, observed value +/- between-participant SE against each model's 5th-95th percentile:

| effect | observed | Model 1 | Model 2 |
| --- | --- | --- | --- |
| A.2 duration benefit on accuracy, at 12% vs the extremes | +0.075 +/- 0.014 | +0.020 to +0.205 | +0.018 to +0.202 |
| duration effect on RT | +49 +/- 6 ms | +2 to +88 | +2 to +90 |
| coherence effect on RT | +33 +/- 6 ms | +2 to +87 | +2 to +79 |
| **B.2 duration effect on the stimulus-period rate** | **+0.057 +/- 0.015 /s** | **-0.241 to -0.004** | -0.125 to +0.143 |
| duration effect on PMET | +49 +/- 6 ms | +4 to +242 | +5 to +768 |
| duration effect on MET | -0.4 +/- 1.4 ms | -183 to -0 | -730 to -1 |
| **C.3 coherence x duration on the delay rate** | **+0.107 +/- 0.051 /s** | **-0.509 to +0.055** | **-0.467 to +0.081** |
| **C.2 delay-period log-rate slope** | **+0.892 +/- 0.092 /s** | **-0.268 to +0.092** | -1.234 to +2.075 |

Three readings worth keeping. **The stimulus-period duration effect excludes Model 1 on its
own**, by sign: its whole admissible range is negative where the observation is positive.
This is Proposition B.2 as a measurement, and it is independent of everything in the delay
period. **The coherence x duration interaction on the delay rate is missed by both models**,
across their whole admissible spaces -- the paper's acknowledged disagreement, now
established as a property of the architectures rather than of a parameter choice. And the
**PMET/MET dissociation is confirmed in the data but is a weak test**: both architectures
can also produce motor-component effects one to two orders of magnitude larger than
anything observed, so finding a small one does not discriminate.

A note on which observed slope to compare against. The mean of the per-cell,
per-participant slopes is 0.892 /s; the slope of the profile pooled across participants is
0.817. The first is the right comparison for a single simulated homogeneous process, and
`frontier.py` uses it: pooling mixes processes at different distances from the gate, and
the log of a mixture is flatter than the mixture of the logs, so the pooled figure
understates any one participant's slope and would make the frontier comparison kinder to
Model 1 than it should be.

### `representative.py` -- the point the illustrative figures are drawn at
Nothing is fitted, so the figures that show simulated dynamics need a rule for choosing one
parameter vector. It is the admissible draw closest to the coordinate-wise median of the
sampled parameters of the draws admissible in BOTH models, each coordinate standardised by
its interquartile range and taken in logarithm where it was sampled in logarithm. The rule
touches only the parameters -- never the rate, the slope, or anything the figures are used
to discuss.

Two traps, both hit and both recorded in the module docstring. The coordinate-wise median
of the *ranges* is not usable: it sits at a standardised distance to the gate near 4, where
the model produces no muscle activity at all. And beta must be set separately, to the median
of the values admissible under Model 2, because requiring admissibility in both models
biases beta downward (0.0075 against 0.0183) -- a large ramp is what produces premature
responses, and the screen removes those -- so taking beta from the draw would illustrate
Model 2 at a ramp smaller than is typical of Model 2, for a reason that has nothing to do
with Model 2.

### `figures_space.py`
`fig_effect_sizes.pdf` is the panel described under `effect_panel.py`.
`fig_frontier.pdf` puts the two quantities that must be reproduced together on the two
axes, draws the analytic frontier as the boundary of what Model 1 can do, and shows where
the sampled draws and the observation fall. `fig_delay_profile.pdf` shows the delay-period
rate curve itself, observed against each architecture, restricted to draws whose mean rate
brackets the observed one -- without that restriction the levels differ by two orders of
magnitude and a bin-by-bin percentile envelope mixes profiles that are not comparable.


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
The figures that show simulated dynamics, all drawn at the representative admissible draw
(`representative.py`): the manuscript's Figure 1 simulated for each model -- the same
parameter vector in both, differing only in beta -- a single trial with the detected
intervals, the rate profiles against data, and the data-only rate/slope relation of
Figure S3. Run `representative.py` first.

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
| Section 3, predictions | derived by hand; checked by `theory_checks.py` |
| Section 3, exclusion of Model 1 | `frontier.py` (analytic), `fused_sweep.py` + `ramp_check.py` (sampled) |
| Section 4, from the process to the observable | `sim.py` (`_bursts`), `theory_checks.py` |
| Section 5, the data | `empirical.py` |
| Sections 6-7, statistical re-analysis | `R/run_stats.R`, `R/emit_tex.R` |
| Section 8, behaviour does not discriminate | derived: Remark A.0 -- beta contributes identically in all six cells, so no behavioural measure can separate the models |
| Section 9, evaluation over the parameter space | `fused_sweep.py`, `effect_sizes.py`, `ramp_check.py` |
| Figures | `figures.py`, `figures_space.py` |
| Appendices A, B, D, E (predictions, no closed form) | derived by hand, not separately checked numerically |
| Appendix C, the rate and the frontier | `theory_checks.py`, `kernel_moments.py`, `frontier.py` |

---

# Things that will bite

**Results are tracked, not scratch.** `results/*.json` (the frontier, the effect-size
summary, the noise decomposition, the representative vector) are the citable numbers and
are committed to git; `results/logs/` and the large `results/*.jsonl` dumps are regenerable
and gitignored (see `.gitignore`), except `ramp_check.jsonl`, which is small and is the
record of the high-precision re-simulation. This
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
simulation and sweep scripts (`sim.py`, `targets.py`, `fused_sweep.py`) extract everything directly from the raw segmented-event and
behavioural files via `empirical.py`/`targets.py` and never touch the JASP exports, so this
bug could not have reached the Monte-Carlo-sweep numbers -- only a
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
follow `sim.py` automatically. Everything downstream of `sim.py`'s `DT`
(`fused_sweep.py`, `export_for_R.py`, `figures.py`) was re-run under the new rate.

**The lateralised readiness potential during the stimulus** reproduces only to about 3%,
in the same direction on all three values. Every effect there is null before and after.

**Nothing here is fitted, and the removal was deliberate.** Two fitting pipelines existed
until 2026-08-31 -- a behaviour-only RT-quantile G2 fit and an all-data extension of it --
and both were removed along with their targets, checkpoints, fit-quality figures and the
G2 machinery. They were only ever diagnostic ("can this architecture reach effects of about
the observed size somewhere in its parameter space"), and the parameter-space evaluation
answers that question over the whole space instead of at one point, without an objective
whose scale is set by a weighting convention and whose model comparison is guaranteed in
advance by nesting. What replaced the behavioural comparison is stronger than it was:
Remark A.0 of the manuscript proves that beta contributes identically in all six cells, so
no behavioural measure in this design can separate the two models -- a theorem where the G2
gap was an observation. If a fit is ever wanted again, it is in the git history, not in the
working tree.

**A statistic that separated the models in theory did not survive contact with the sweep,
and is gone.** The curvature of the delay-period profile was added and removed the same
day. The theory is sound -- Model 1's log-rate slope decays geometrically at 2*kappa,
Model 2's linearly, and the ratio of the first-quarter to the last-quarter slope separates
them 10-18x against 1.2-3x -- but a fitted quadratic coefficient measures something else:
it equals -k^2/2 with k = alpha/sigma_V, so it confounds the SHAPE of the ramp with its
SIZE, and a model reproducing a large slope necessarily produces a large curvature. Worse,
there was no matched comparison available, because Model 1 never reaches the observed slope
at any rate. Two further traps are worth recording because they are easy to fall into
again: pooling across participants FLATTENS curvature (log-sum-exp of concave functions is
not concave -- measured, -0.095 per participant against -0.034 pooled, a factor of three in
the direction one might assume impossible); and the last quarter of the window is depressed
in the model by an effect that is not curvature at all, since `measures()` and
`empirical.py` both score the burst preceding the response as the response rather than as a
partial activation, so a trial whose final approach begins before the reaction signal loses
an event from the end.
