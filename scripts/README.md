# Analysis code for the theoretical note

Two languages, split by job.

**Python** owns the model: the stochastic simulation, the electromyographic detector, the
parameter-space evaluation, the closed-form rate formulas, the extraction of events from
the archive, and the figures. Nothing is fitted -- the models' predictions are derived
analytically and their effect-size ranges come from sampling the parameter space, not from
estimating a parameter vector.

**R** owns the statistics that this project computes: the analyses of variance on the
electromyographic and electroencephalographic measures and the trend decompositions -- the LEVEL, TREND and
QUAD contrasts on the binned measures, which since 2026-09-07 are the only statistics the
paper reports on the partial-activation counts; the unbinned point-process analysis was
removed with the cold start (see below).

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
  It had never fired because `--with-sim` had never run to completion. (That flag, and the
  analysis of simulated counts it switched on, were removed on 2026-09-03; see
  `R/run_stats.R`.)

The two sides meet at CSV files written by `export_for_R.py`, which carry counts and
exposures for the data and for each model simulated at the representative parameter
vector, so that the identical statistical analysis runs on the empirical side.

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
Rscript run_stats.R
```

On a machine that already has `renv.lock` (checked in after the first `renv::snapshot()`),
skip straight to `Rscript -e 'renv::restore()'` to reproduce the exact recorded package
versions instead of re-resolving them.

## Order of execution

All paths below are relative to `scripts/` (run `cd scripts` first); `results/` and
`plots/` already exist in the repository. Prefix any of these with `uv run python` to use
the managed virtualenv. The pipeline runs top to bottom: everything downstream of the sweep
reads its dump.

Everything downstream of the sweep is wrapped in `./rerun.sh`, which runs the steps below
in order and takes about two minutes. The sweep is not in it: it is long, and is run
sharded across cores by `./run_sweep.sh`.

**The stationary start (2026-09-07).** `sim.py` now draws the maintenance noise at
maintenance onset from its own stationary law rather than setting it to zero
(`WARM_START`). One line, and it changes what the delay-period analysis can claim. Under
the cold start the fluctuation of the signal about the maintained level carried a mode
decaying at 2*kappa -- as slow as 0.4 /s over the sampled range -- so the transient
outlasted the delay and every delay-period statement was conditional on it. Started at
stationarity the transient decays at a rate of at least lambda >= 10 /s, so the analysis
window (which opens 1 s after stimulus offset) lies wholly in the stationary regime, and
Model 1's delay-period rate is exactly constant at every parameter setting. Two
consequences run through this directory: the frontier has no content left, and the
discriminating statistic is the raw-rate TREND rather than the log-rate slope. Cold-start
dumps and outputs are kept under their unsuffixed names for comparison.

```
python delaymodel/empirical.py    -> (printout)              extraction, checked against the archive
python delaymodel/targets.py      -> (printout)              empirical summaries per cell

# 1. the parameter-space evaluation  (~12 min over six shards; see fused_sweep.py below)
#    ./run_sweep.sh does this natively, one shard per core, resumable
SHARD_N=6 SHARD_I=$i \
  python fused_sweep.py 1200 4000 results/shards/shard_$i.jsonl
python sweep_summary.py            -> results/sweep_screened.jsonl  the concatenated dump
                                     results/sweep_summary.json

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

# 2b. leak robustness (reviewer C2.3): the same draws at a non-zero leak
for G in 0.25 0.50 1.00; do tag=${G/./}
  for i in 0 1 2 3; do
    GAMMA=$G OUT=results/leak/leak_g${tag}_sh$i.jsonl SHARD_N=4 SHARD_I=$i \
      python leak_sweep.py > results/logs/leak/g${tag}_sh$i.log 2>&1 &
  done; wait
  cat results/leak/leak_g${tag}_sh*.jsonl > results/leak/leak_g${tag}.jsonl
done
python leak_summary.py            -> results/leak_summary.json     the table quoted in the letter

# 3. the analytic checks, which need no sampling at all
python mkfig.py                   -> plots/fig_approx_error.pdf    error of the stationary
                                                                   approximation over the
                                                                   sampled domain (manuscript
                                                                   Figure 8)
python theory_checks.py           -> (printout)                    verification of the Appendix

# 4. the CSV bridge, then the figures
# order matters: figures.py reads results/representative.json, figures_space.py reads
# pma_counts.csv and results/effect_panel.json
python representative.py          -> results/representative.json   the vector Figure 2B is
                                                                   drawn at
python export_for_R.py            -> pma_counts.csv                counts and exposures, empirical
                                                                   (written at scripts/ root: R's
                                                                   SIMDIR default, "..", expects it
                                                                   there)
python figures.py                 -> plots/fig_traj_m{2,1}.pdf     manuscript Figure 1
                                     plots/fig_single.pdf         manuscript Figure 2B
python effect_panel.py            -> results/effect_panel.json    effect-size ranges vs observed
python figures_space.py           -> plots/fig_effect_sizes.pdf   those ranges, as a panel
                                     plots/fig_rate_trend.pdf     rate against trend, both models
                                     plots/fig_delay_profiles.pdf the delay rate curve

# 6. the manuscript's own copies, at 200 dpi (there is no script for this; the mapping is
#    fig01,02 <- fig_traj_m2,m1;  fig04 <- fig_single;  fig08 <- fig_rate_trend;
#    fig09 <- fig_delay_profiles;  fig10 <- fig_effect_sizes;  fig11 <- fig_approx_error.
#    fig03, 05, 06 and 07 are not produced here.)
cd plots && pdftoppm -r 200 -png -singlefile fig_rate_trend.pdf ../../manuscript/latex/figures/fig08

# 5. the statistics, in R
Rscript R/run_stats.R             -> R/stats_tables.txt            all tables (empirical)
                                     R/stats_results.csv           one row per effect
                                                                   (includes the behavioural
                                                                   ANOVAs of Table S1
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
and stimulus partial-activation rates, the delay-period TREND and QUAD contrasts, the
delay log-rate slope, premature-response rate. These are descriptive targets used by the
figures and by the screen's observed-rate window; they are not fitted to.

TREND and QUAD come from the archive's own bin export (`analysis_PMA/data/
PMA_rate_delay_JASP.csv`) rather than from the event-file reconstruction that produces
everything else here. That export is the file the published analysis ran on, so contrasts
taken from it reproduce the manuscript's F values exactly; the reconstruction gives a mean
trend of 0.5551 against the export's 0.5604, and a 1% offset between the comparison target
and the number in the Results is not worth having.

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

Latin-hypercube sample over the ranges of Table 1; the production run uses 1200 draws at
4000 trials per cell per model per draw (~6.5 s/draw, ~25 min over six shards). At every draw:
the admissibility screen, the two RT-effect criteria (coherence, duration), the
delay-period criteria on the raw-rate TREND and on the log-rate slope (all six cells
significantly positive, `trend_criterion` and `slope_criterion`) and the
motor-execution-time invariance criterion. The same sampled vector is evaluated in both
models, with beta set to zero in Model 1, so the two are compared at identical parameters.

Resumable by line count of `OUT`. Shardable: `SHARD_N`/`SHARD_I` in the environment make
shard j take the draws i with i % n == j, each writing its own file. The LHS cube is
regenerated identically from the fixed seed in every process, so `cat`-ing the shard files
reproduces the unsharded run exactly -- verified by re-running the sweep sharded and
recovering identical per-draw statistics.

```
for i in 0 1 2 3 4 5; do
  SHARD_N=6 SHARD_I=$i uv run python3 fused_sweep.py 1200 4000 \
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
    delay) to 33 ms. The range is wide because nothing constrains kappa, and since the
    stationary start it favours neither model: the transient decays at at least lambda
    whatever kappa is, so no value of kappa gives Model 1 a delay-period ramp (Appendix,
    Corollary A.2). Under the cold start the same width was load-bearing in Model 1's
    favour -- small kappa was its only route to a rate still rising at the end of the
    delay -- and a floor at kappa = 2.0 imposed on 2026-08-30 was reverted the next day for
    that reason. The floor's justification, that the stationary crossing-rate approximation
    loses validity below it, is now simply true and harmless: the approximation is used,
    and `mkfig.py` measures its error over the whole domain (1.4e-5 inside the window).
  - `lam` in [10.0, 50.0], log-uniform. 1/lam is the time constant of the filter relating
    decision-related input to motor preparation, so the range corresponds to 100 down to
    20 ms; lam = 100 would mean 10 ms, faster than plausible neuromuscular filtering. This
    bound is a modelling commitment rather than a convenience, and is stated as one in the
    manuscript's Table 1 note: lam sets how fast the fluctuation settles, so it is what
    puts the analysis window in the stationary regime (Appendix, Remark A.1), and a range
    chosen after seeing that consequence would be circular. It was fixed from the
    neuromuscular argument before any of this appendix existed.

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
reached the observed MAGNITUDE of the stimulus-period duration effect. That comparison
therefore has to be made on the SIGN, which is what the proposition predicts anyway -- and
there the separation is undamaged, 33% against 0%.

**Result** (`results/sweep_screened.jsonl`, 1200 draws, NTR 4000). Model 1: 160 admissible;
delay-period trend median -0.007 /s^2, 5th-95th percentile -0.058 to +0.004, largest value
+0.007; **0** draws with all six condition trends significantly positive and **0** anywhere
in the observed 95% confidence interval [0.375, 0.745]. Model 2: 145 admissible; trend
median +0.052, -0.017 to +0.418, largest +1.157; 78 with all six positive and 6 in the
observed interval with the rate also in range. The observed +0.5604 sits at Model 2's
96.6th percentile and above every draw of Model 1.

The fraction of admissible settings carrying the OBSERVED SIGN is the compact summary of the
whole comparison, and `effect_panel.py` reports it:

| effect | Model 1 | Model 2 |
| --- | --- | --- |
| the seven effects both architectures share | 86-98% | 88-94% |
| **duration effect on the stimulus-period rate** | **0%** | **33%** |
| **the delay-period rate trend** | **1%** | **80%** |
| coherence x duration on the delay rate | 3% | 5% |

Two effects separate the architectures, one is missed by both, seven are indistinguishable.

**Proposition B.2 was half wrong and has been corrected.** It claimed that Model 2 predicts a
higher stimulus-locked rate after the SHORT stimulus and Model 1 the opposite. Model 1's half
holds at 96% of its admissible settings. Model 2's does not: temporal preparation opposes the
level effect but overturns it at only 33%, because the ramp's size and the distance to the
gate are independent parameters. The test is therefore ONE-SIDED -- it refutes Model 1 and
does not confirm Model 2 -- and the proposition, the Results paragraph and Remark B.3 now say
so. This was caught by reading the effect panel, which showed both models' medians on the
same side of zero where the text said they should be on opposite sides.

A fresh run with a different seed will not reproduce the exact admissibility counts; the
criteria and their signs are robust, the counts are not guaranteed identical run to run.
Keep the summary JSON alongside the script rather than regenerating it casually.

### Retired: `kernel_moments.py` and `frontier.py` (2026-09-07)

Both computed properties of the COLD-START model, in which the maintenance noise began at
zero and the fluctuation of the signal about the maintained level was still filling in
during the delay window. `kernel_moments.py` held the exact non-stationary moments and the
Kac-Rice rate that follows from them; `frontier.py` used them to compute the ramp-rate
frontier -- the largest delay-period log-rate slope Model 1 could reach at each level of
the rate -- which was the analytic core of the delay-period argument and the manuscript's
old Figure 5.

The stationary start abolishes the trade-off they described. There is no filling-in inside
the analysis window, Model 1's delay-period rate is exactly constant (Appendix, Corollary
A.2), the largest trend the formulas return anywhere in the sampled ranges is 3e-11 /s^2,
and a frontier of that is a flat line at zero. The manuscript's Figure 5 is now the
rate-against-trend scatter (`figures_space.fig_rate_trend`), which needs no envelope
because the boundary of what Model 1 can do is the horizontal axis. The exact moments that
`kernel_moments.py` held survive where they are still needed, in `mkfig.py`, which uses
them for one purpose only: to measure how far the stationary approximation is from the
exact rate over the sampled domain (1.4e-5 inside the window, manuscript Figure 8).

Both files are deleted rather than kept as dead code; they are in the git history at
commit `1e361ab` if the cold-start argument ever needs to be reread.

### `ramp_check.py` -- the existence question, at high precision
An existence question is badly served by raising the trial count uniformly across a sweep:
the statistic that answers it is a MAXIMUM over draws, and a maximum over noisy estimates
is biased upward, so precision belongs on the candidates rather than spread evenly. Stage 1
is `fused_sweep.py` at moderate NTR, to find candidates; this is stage 2, re-simulating the
best of them at NTR 40000 on several independent seeds, selected two ways -- the best among
those producing a realistic rate (the honest conditional question) and the best overall
regardless of rate.

Since the stationary start the question has an exact answer -- Corollary A.2: Model 1's
delay-period rate is constant, so its trend is zero everywhere -- and the script's job is
to check that the simulation agrees once the Monte-Carlo error is small enough to see the
difference. It reports the verdict on exit, recomputed from the whole file so an
interrupted run still says something.

Result (`results/ramp_check.jsonl`), seven candidate draws at NTR 40000, three seeds each:

| draw | rate /s | mean trend /s^2 | SE |
|------|---------|-----------------|-------|
| 628  | 0.216   | +0.0007         | 0.0005 |
| 461  | 0.217   | +0.0000         | 0.0005 |
| 6    | 0.141   | -0.0002         | 0.0006 |
| 1164 | 0.120   | -0.0014         | 0.0009 |
| 15   | 0.972   | -0.0022         | 0.0002 |
| 220  | 0.216   | -0.0034         | 0.0007 |
| 965  | 1.798   | -0.0085         | 0.0028 |

The largest is +0.0007 /s^2 against an observed +0.5604: eight hundred times below, and
zero to within the error. The negative values are the response-burst exclusion of Appendix
Remark A.5, and they scale with the rate exactly as that account requires -- the most
active draw is the most negative. This is the independent confirmation that the sweep's
+0.0074 maximum is noise and not a residual ramp.

### `leak_sweep.py` -- leak robustness (reviewer C2.3)

Re-runs, at a non-zero leak on the decision variable, exactly the draws that were admissible
in either model at zero leak. It imports `fused_sweep` and reuses its cube, `to_param`,
`run_model` and `admissible`, so the parameter vectors are identical to the main evaluation
and only `gamma` differs; `results/leak/adm_idx.json` is the list of those 220 draw indices,
written by `leak_index.py` from `results/sweep_screened.jsonl`. Resumable by line count of
`OUT`; `SHARD_N`/`SHARD_I` split the draw list across cores. Records have the same shape as
`sweep_screened.jsonl` minus the `raw` block, plus a criterion block per model for the
TREND and for the log-rate slope. `run_leak.sh` drives all three leak rates one shard per
core; `chunk_leak.sh` does the same in bounded slices, for a shell that gets a fresh
process and a hard time limit per call. `leak_summary.py` reduces them to the table.

Result (`results/leak_summary.json`), on the 160 draws admissible for Model 1 and the 145
for Model 2 at gamma = 0:

| gamma | M1 adm. | M2 adm. | RT coherence | RT duration | M2 all six delay trends > 0 | M1 all six > 0 |
|-------|---------|---------|--------------|-------------|-----------------------------|----------------|
| 0.00  | 160     | 145     | 136/145      | 136/145     | 78                          | 0              |
| 0.25  | 155     | 143     | 137/145      | 137/145     | 76                          | 0              |
| 0.50  | 152     | 142     | 136/145      | 137/145     | 71                          | 0              |
| 1.00  | 138     | 130     | 138/145      | 137/145     | 77                          | 0              |

**Model 1 never produces all six condition trends significantly positive, at any leak
rate.** That is worth more than the same sentence at zero leak: Corollary A.2 proves the
delay-period rate constant for the model as specified, and a leak on the decision variable
is a perturbation the corollary does not cover, so this is an independent check rather than
a restatement. Model 2 clears the criterion at roughly half its admissible draws
throughout. The log-rate slope column, recorded alongside, agrees to within three draws
everywhere.

Nothing directional moves. What erodes is admissibility itself, and the draws that drop out
are those with the smallest `k*c/theta` (median 5.4 against 7.8 for those that survive): a
leak caps the accumulator at `k*c/gamma`, and once that asymptote approaches `theta` the
threshold stops being reached. That ratio is the criterion for when a leak can matter at
all; its smallest value at 50% coherence over the admissible draws is 1.5 s^-1, above the
largest rate tested. The numbers are quoted in the response letter (C2.3), not in the
manuscript.

(The cold-start version of this table -- 164/148 admissible, the criterion on the log-rate
slope -- is in the git history at commit `1e361ab` and earlier. It said the same thing.)

### `effect_sizes.py` -- the citable summary, with its noise measured
Stage 1 distils the sweep dump into the counts and distributions quoted in the manuscript.
Stage 2 measures how much of the reported spread is simulation noise rather than parameter
sensitivity, by re-simulating a subset of admissible draws on several seeds and decomposing
the variance of the per-draw trend into between-draw and within-draw (Monte-Carlo) parts.

Since 2026-09-07 the statistic throughout is TREND -- the linear orthogonal polynomial
contrast of the raw rate, in /s^2 -- and not the offset-Poisson log-rate slope. The
log-rate slope is a ratio, so a draw can earn a steep one by ramping from almost nothing;
the trend is the product of the level and the relative rise, so it demands both, which is
what the observation shows. The acceptance band is the between-participant 95% confidence
interval on the observed trend, [0.375, 0.745] /s^2, and not the spread of the six
condition means ([0.469, 0.634], still reported alongside): a spread across conditions is
not an uncertainty on the quantity a model is asked to reproduce. Both counts are in the
output, as are two prior-free summaries -- `pct_of_observed` (Model 2: 96.6) and
`ratio_observed_to_max` (Model 1: 75) -- which are what the Results lead with, because a
count of draws in a band reports the Latin-hypercube prior as much as the architecture.
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

Result, observed value +/- between-participant SE against each model's 5th-95th percentile
(`results/effect_panel.json`; the panel letters are the manuscript's Figure 10 A-I):

| effect | observed | Model 1 | Model 2 |
| --- | --- | --- | --- |
| A duration benefit on accuracy, at 12% vs the extremes | +0.075 +/- 0.014 | +0.015 to +0.196 | +0.017 to +0.204 |
| B duration effect on RT | +49 +/- 6 ms | +1 to +80 | +1 to +85 |
| C coherence x duration on RT | +15 +/- 8 ms | +1 to +78 | +1 to +77 |
| D coherence effect on the stimulus-period rate | +0.054 +/- 0.020 /s | +0.012 to +0.686 | +0.001 to +0.416 |
| **E duration effect on the stimulus-period rate** | **+0.057 +/- 0.015 /s** | **-0.340 to -0.007** | -0.217 to +0.043 |
| F duration effect on PMET | +49 +/- 6 ms | +2 to +108 | +2 to +129 |
| G duration effect on MET | -0.5 +/- 0.8 ms | -46 to -0.2 | -60 to -0.7 |
| **H coherence x duration on the delay rate** | **+0.107 +/- 0.051 /s** | **-0.769 to +0.007** | **-0.750 to +0.008** |
| **I delay-period rate trend** | **+0.560 +/- 0.094 /s^2** | **-0.058 to +0.004** | -0.017 to +0.418 |

Three readings worth keeping. **The stimulus-period duration effect excludes Model 1 on its
own**, by sign: its whole admissible range is negative where the observation is positive,
and it is independent of everything in the delay period. **The coherence x duration
interaction on the delay rate is missed by both models**, across their whole admissible
spaces -- the paper's acknowledged disagreement, established as a property of the
architectures rather than of a parameter choice. And the **PMET/MET dissociation is
confirmed in the data but is a weak test**: both architectures can also produce
motor-component effects far larger than anything observed, so finding a small one does not
discriminate.

A note on the observed value to compare against, which the switch to the TREND removed a
trap from. The log-rate slope had two of them -- 0.892 /s as the mean of the per-cell,
per-participant slopes and 0.817 as the slope of the pooled profile -- because the log of a
mixture is flatter than the mixture of the logs, so which one to use had to be argued. The
trend is linear in the rate, so the two coincide exactly at +0.5604 /s^2 and there is
nothing to argue.

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
`fig_effect_sizes.pdf` is the panel described under `effect_panel.py`. `fig_rate_trend.pdf`
puts the two quantities that must be reproduced together on the two axes -- the delay-period
rate and its linear trend -- with the admissible draws of both models and the observation;
it needs no computed frontier, because Model 1's prediction is the horizontal axis.
`fig_delay_profiles.pdf` shows the delay-period rate curve itself, each draw divided by its
own mean over the window so that draws whose rates differ by two orders of magnitude are
comparable on shape.


### `export_for_R.py`
Writes the bridge file. `pma_counts.csv` has one row per participant, cell, period and
bin, with the raw count and the number of trials contributing to that bin. Until
2026-09-05 it also wrote `sim_counts_m{1,2}.csv`, the same quantities from each model at
the representative vector; nothing read them -- R analyses the empirical counts only, and
what the models predict is reported as effect-size ranges over the parameter space -- so
the files and the function that built them were removed.

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
The figures that show simulated dynamics: the manuscript's Figure 1, drawn for each model
over thirty admissible draws -- the same thirty in both, differing only in beta -- and
Figure 2B, one simulated trial with the detected intervals, at the representative vector.
Run `representative.py` first (Figure 1 does not need it; Figure 2B does).

Two functions were removed on 2026-09-07. `fig_rates` drew the model's rate profile against
the data at that one vector, and is superseded by `figures_space.fig_delay_profiles`, which
shows the whole admissible set. `fig_shape` drew the log-rate slope against the rate across
the six conditions; its manuscript figure went when the log-rate slope did.

### `theory_checks.py`
Verification of the Appendix (the Kac-Rice rate derivation; renumbered from Appendix A in
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
(trend contrasts). Writes `stats_tables.txt` and `stats_results.csv`. It analyses the
empirical counts only; the `--with-sim` flag and the analysis of simulated counts were
removed on 2026-09-03, because a parameter vector is one hypothetical participant and the
simulated design therefore has no between-participant factor.

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
| Section 3, exclusion of Model 1 | derived: Corollary A.2 (analytic); `fused_sweep.py` + `ramp_check.py` (sampled) |
| Section 4, from the process to the observable | `sim.py` (`_bursts`), `theory_checks.py` |
| Section 5, the data | `empirical.py` |
| Sections 6-7, statistical re-analysis | `R/run_stats.R`, `R/emit_tex.R` |
| Section 8, behaviour does not discriminate | derived: Remark A.0 -- beta contributes identically in all six cells, so no behavioural measure can separate the models |
| Section 9, evaluation over the parameter space | `fused_sweep.py`, `effect_sizes.py`, `ramp_check.py` |
| Figures | `figures.py`, `figures_space.py` |
| Appendices A, B, D, E (predictions, no closed form) | derived by hand, not separately checked numerically |
| The Appendix, the rate and the stationary approximation | `theory_checks.py`, `mkfig.py` |

---

# Things that will bite

**Results are tracked, not scratch.** `results/*.json` (the effect-size summary, the effect
panel, the noise decomposition, the leak summary, the representative vector) are the
citable numbers and are committed to git; `results/logs/`, `results/shards/`, `plots/` and
the large `results/*.jsonl` dumps are regenerable and gitignored (see `.gitignore`), except
`ramp_check.jsonl`, which is small and is the record of the high-precision
re-simulation. This
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
fresh, un-cached run of the R trend-contrast reanalysis (Tables S1-S4) was at
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

**The curvature is an exact prediction that the simulation cannot measure cleanly, and the
manuscript says so.** Under the stationary start the two models make opposite and exact
claims about the quadratic contrast of the delay-period rate: zero under Model 1, since a
constant rate has zero curvature, and positive under Model 2 wherever the gate is more than
one fluctuation standard deviation away, since d2(nu)/drho2 = nu (beta/delta sigma_H)^2
(a^2 - 1). The data give +0.2099 /s^3, F(1,21) = 18.891. The simulated curvatures do not
reach it: 0 of 160 Model 1 draws and 19 of 145 Model 2 draws have all six condition
curvatures significantly positive, and both medians sit near zero.

The reason is measurement, not model, and it is worth recording because it is easy to
misread as a failed prediction. `measures()` and `empirical.py` both score the burst
preceding the response as the response rather than as a partial activation, so a trial
whose final approach begins before the response signal loses an event from the END of the
window -- and the quadratic contrast weights the ends. Restricted to the first eighteen
bins Model 2's median curvature turns positive and Model 1's halves. The Results therefore
report the curvature as the analytic prediction it is and the simulated distribution as a
lower bound on it.

An older trap, from the log-rate era, is still worth knowing: pooling across participants
FLATTENS curvature (log-sum-exp of concave functions is not concave -- measured, -0.095 per
participant against -0.034 pooled, a factor of three in the direction one might assume
impossible). The raw-rate contrasts do not have this problem, being linear in the rate.
