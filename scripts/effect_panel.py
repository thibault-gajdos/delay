"""
effect_panel.py -- the effect-size ranges each architecture can produce, on every
prediction the paper tests, against the observed value.

This is what the parameter-space evaluation is FOR. The predictions of Appendices A-E are
ordinal and hold at every admissible setting; what a derivation cannot supply is a sense of
scale. For each effect the paper tests, this reports the distribution of the effect over
the admissible draws of each model -- median and 5th-95th percentile -- next to the value
observed in the data. A model whose range excludes the observation is quantitatively
inadequate on that effect, whatever the sign of its prediction.

The effects, all computed from the same dump and in the units the Results report them in:

  accuracy, 50% coherence      mean over the two 50%-coherence cells
  coherence effect on RT       mean RT at 0% minus at 50% coherence, marginal over duration
  duration effect on RT        mean RT after the 0.3 s minus after the 2 s stimulus
  MET range across cells       max minus min of cell-mean motor execution time -- a
                               prediction of INVARIANCE, so a small range is the prediction
                               and a large one would falsify both models
  delay-period PMA rate        mean over the 1-3 s window
  delay-period log-rate slope  offset-Poisson slope over the same window

Writes results/effect_panel.json.  Usage: python3 effect_panel.py [DUMP.jsonl]
"""
import json, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.targets import build

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
OUT = 'results/effect_panel.json'
# cell order (targets.CELLS): (c, D) = (0,.3) (.12,.3) (.5,.3) (0,2) (.12,2) (.5,2)
SHORT, LONG = [0, 1, 2], [3, 4, 5]
C00, C50 = [0, 3], [2, 5]


def effects(rt, acc, met, rate, slope):
    """the six effects, from per-cell arrays"""
    return dict(
        acc50=float(np.nanmean(acc[C50])),
        coh_rt=float(np.nanmean(rt[C00]) - np.nanmean(rt[C50])) * 1000.0,
        dur_rt=float(np.nanmean(rt[SHORT]) - np.nanmean(rt[LONG])) * 1000.0,
        met_range=float(np.nanmax(met) - np.nanmin(met)) * 1000.0,
        rate=float(np.nanmean(rate)),
        slope=float(np.nanmean(slope)))


NAMES = dict(acc50='accuracy at 50% coherence',
             coh_rt='coherence effect on RT (ms)',
             dur_rt='duration effect on RT (ms)',
             met_range='MET range across cells (ms)',
             rate='delay-period PMA rate (/s)',
             slope='delay-period log-rate slope (/s)')

recs = [json.loads(l) for l in open(DUMP)]

E = build()
# the observed delay-period rate per cell, in the same units the simulation reports:
# counts per 0.1 s bin summed over the 20 bins of the 2 s window, divided by 2 s
obs_rate_cells = E['dly'].mean(0).sum(axis=-1) / 2.0
obs = effects(E['rt'].mean(0), E['acc'].mean(0), E['met'].mean(0),
              obs_rate_cells, np.nanmean(E['slope'], 0))
# between-participant SE of each effect, so the reader can see what "the observed value"
# is worth before comparing it with a model range
per_sub = [effects(E['rt'][s], E['acc'][s], E['met'][s],
                   E['dly'][s].sum(axis=-1) / 2.0, E['slope'][s])
           for s in range(E['rt'].shape[0])]
obs_se = {k: float(np.nanstd([p[k] for p in per_sub], ddof=1) / np.sqrt(len(per_sub)))
          for k in NAMES}

out = dict(source=DUMP, names=NAMES, observed=obs, observed_se=obs_se, models={})
for tag in ['model1', 'model2']:
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    rows = [effects(np.asarray(m['raw']['rt_mean'], float), np.asarray(m['raw']['acc'], float),
                    np.asarray(m['raw']['met_mean'], float), np.asarray(m['raw']['rate'], float),
                    np.asarray(m['raw']['slope'], float)) for m in adm]
    out['models'][tag] = dict(n=len(rows))
    for k in NAMES:
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
        out['models'][tag][k] = dict(
            median=float(np.median(v)), p5=float(np.percentile(v, 5)),
            p95=float(np.percentile(v, 95)), min=float(v.min()), max=float(v.max()),
            covers_observed=bool(np.percentile(v, 5) <= obs[k] <= np.percentile(v, 95)))

json.dump(out, open(OUT, 'w'), indent=1)
w = max(len(v) for v in NAMES.values())
print(f"{'effect':{w}} {'observed':>12} {'Model 1  5-95%':>26} {'Model 2  5-95%':>26}")
for k, name in NAMES.items():
    o = obs[k]; se = obs_se[k]
    f = '%.3f' if k in ('acc50', 'rate', 'slope') else '%.0f'
    cells = []
    for tag in ['model1', 'model2']:
        r = out['models'][tag][k]
        mark = ' ' if r['covers_observed'] else '*'
        cells.append(f"{f % r['p5']} to {f % r['p95']}{mark}")
    print(f"{name:{w}} {f % o} +/-{f % se}  {cells[0]:>25} {cells[1]:>25}")
print("\n* = the model's 5th-95th percentile does not cover the observed value")
