"""
effect_sizes.py -- the citable effect-size summary for the screened sweep, with the
Monte-Carlo component of its variability measured rather than assumed.

Two stages, deliberately separated:

  STAGE 1 (summary).  Reads the sweep dump and reports, per model, how many draws pass the
  screen and what the delay-period TREND does across them -- the linear orthogonal
  polynomial contrast of the raw rate against time, which is the statistic the data
  analysis uses (Methods: TREND).  Median, 5-95% range, maximum, and how many draws land
  inside the observed per-condition band [OBS_LO, OBS_HI] = [0.469, 0.634] /s^2 while
  ALSO producing a delay rate inside the range
  real participants span.  The joint requirement is the point: a draw that ramps at
  1/300th of the observed activation rate has not reproduced the observation.

  STAGE 2 (noise).  A single draw's trend is itself a noisy estimate, so the spread
  reported in stage 1 mixes real parameter sensitivity with simulation noise.  A subset of
  admissible draws is therefore re-simulated at N_SEEDS independent seeds and the variance
  of the per-draw mean trend decomposed:  total = between-draw + within-draw (MC).  The
  within-draw component is what an increase in NTR would remove; the between-draw
  component is the effect size the sweep exists to measure.  Reported as SDs, with the MC
  share of total variance.

  Note on maxima.  Stage 1's `max` is a maximum over noisy draws and is therefore biased
  upward -- the winner's-curse point.  Under the stationary start this matters especially
  for Model 1, whose delay-period trend is ZERO analytically (largest attainable value over
  the sampled ranges: 3e-11 /s^2), so its simulated maximum is pure Monte-Carlo noise.  It is reported as an upper bound on what the model
  can do, never as an estimate of what it does, and any draw that clears a threshold on
  this statistic must be re-simulated at high NTR (ramp_check.py) before it is believed.

Usage: python3 effect_sizes.py DUMP.jsonl [N_SEEDS] [NTR] [N_SUBSET]
"""
import json, sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened_warm.jsonl'
N_SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
NTR = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
N_SUBSET = int(sys.argv[4]) if len(sys.argv) > 4 else 40
SHARD_N = int(os.environ.get('ES_SHARD_N', '1'))
SHARD_I = int(os.environ.get('ES_SHARD_I', '0'))

# fused_sweep.py parses sys.argv at import time (it is primarily a CLI script), so hide
# this script's arguments from it rather than have it try to read DUMP as N_DRAWS.
_argv, sys.argv = sys.argv, sys.argv[:1]
import fused_sweep as F
sys.argv = _argv
# Keep the cold-start outputs: a warm dump writes to *_warm.json so the previously
# committed cold-start files are never overwritten.
_SFX = '_warm' if 'warm' in os.path.basename(DUMP) else ''
OUT = (f'results/effect_sizes_stage1{_SFX}.json' if N_SUBSET <= 0 else
       f'results/effect_sizes{_SFX}_sh{SHARD_I}.json' if SHARD_N > 1
       else f'results/effect_sizes{_SFX}.json')

# Observed per-condition delay-period TREND (raw-rate linear contrast, /s^2).  Unlike the
# log-rate slope this statistic is LINEAR in the rate, so the mean of the per-cell,
# per-participant contrasts and the contrast of the pooled profile coincide exactly
# (+0.5604) -- the 0.892-vs-0.817 ambiguity of the log-rate slope does not arise.
OBS_LO, OBS_HI = 0.4685, 0.6342
OBS_TREND = 0.5604
OBS_RATE, OBS_RMIN, OBS_RMAX = 0.709, 0.085, 2.785   # window-mean rate, /s

recs = [json.loads(l) for l in open(DUMP)]
summary = dict(source=DUMP, n_draws=len(recs), ntr_stage1=None,
               obs=dict(trend_lo=OBS_LO, trend_hi=OBS_HI, trend_mean=OBS_TREND,
                        rate_mean=OBS_RATE, rate_min=OBS_RMIN, rate_max=OBS_RMAX))


def q(x, p):
    return float(np.percentile(x, p)) if len(x) else float('nan')


for tag in ['model1', 'model2']:
    adm = [r for r in recs if r[tag]['admissible']]
    sl = np.array([np.nanmean(np.asarray(r[tag]['raw']['trend'], float)) for r in adm])
    rt = np.array([np.nanmean(np.asarray(r[tag]['raw']['rate'], float)) for r in adm])
    ok = np.isfinite(sl)
    sl_f, rt_f = sl[ok], rt[ok]
    inband = (sl_f >= OBS_LO) & (sl_f <= OBS_HI)
    rate_ok = (rt_f >= OBS_RMIN) & (rt_f <= OBS_RMAX)
    row = dict(
        n_admissible=len(adm),
        n_with_measurable_trend=int(ok.sum()),
        trend=dict(median=q(sl_f, 50), p5=q(sl_f, 5), p95=q(sl_f, 95),
                   max=float(sl_f.max()) if len(sl_f) else float('nan')),
        rate=dict(median=q(rt_f, 50), p5=q(rt_f, 5), p95=q(rt_f, 95)),
        n_trend_in_observed_band=int(inband.sum()),
        n_trend_in_band_and_rate_observed=int((inband & rate_ok).sum()),
        # recomputed from the raw-rate contrast rather than read from the sweep record's
        # `slope` block, which holds the Poisson log-rate criterion
        n_all_six_trends_positive=sum(
            1 for r in adm
            if np.all(np.isfinite(np.asarray(r[tag]['raw']['trend'], float)))
            and np.all(np.isfinite(np.asarray(r[tag]['raw']['trend_se'], float)))
            and np.all(np.asarray(r[tag]['raw']['trend'], float)
                       - 1.96 * np.asarray(r[tag]['raw']['trend_se'], float) > 0)),
        n_coh_supported=sum(1 for r in adm if r[tag]['rt'] and r[tag]['rt']['coh_supported']),
        n_dur_supported=sum(1 for r in adm if r[tag]['rt'] and r[tag]['rt']['dur_supported']),
    )
    summary[tag] = row

# ---------------------------------------------------------------- stage 2: MC noise
F.NTR = NTR
rng = np.random.default_rng(20260831)
for tag, model in ([] if N_SUBSET <= 0 else [('model1', 1), ('model2', 2)]):
    adm = [r for r in recs if r[tag]['admissible']]
    if not adm:
        continue
    pick = list(rng.choice(len(adm), size=min(N_SUBSET, len(adm)), replace=False))
    pick = [j for n, j in enumerate(pick) if n % SHARD_N == SHARD_I]
    per_draw = []
    for j in pick:
        r = adm[j]
        vals = []
        for s in range(N_SEEDS):
            res = F.run_model(r['params'], r['alpha'] if model == 2 else 0.0, model,
                              900000 + 1000 * r['draw'] + 137 * s)
            vals.append(float(np.nanmean(np.asarray(res['trend'], float))))
        per_draw.append(vals)
        print(f"{tag} draw {r['draw']}: " +
              " ".join(f"{v:+.3f}" for v in vals), flush=True)
    A = np.array(per_draw)                       # n_draws x n_seeds
    good = np.isfinite(A).all(axis=1)
    A = A[good]
    means = A.mean(axis=1)
    within = A.var(axis=1, ddof=1).mean()        # MC variance per draw
    total = means.var(ddof=1) + within * (1 - 1.0 / N_SEEDS)
    between = max(total - within, 0.0)
    summary[tag]['noise'] = dict(
        n_draws_resimulated=int(A.shape[0]), n_seeds=N_SEEDS, ntr=NTR,
        sd_total=float(np.sqrt(total)), sd_mc=float(np.sqrt(within)),
        sd_between_draw=float(np.sqrt(between)),
        mc_share_of_variance=float(within / total) if total > 0 else float('nan'))

text = json.dumps(summary, indent=1)
print(text)
open(OUT, 'w').write(text)
