"""leak_summary.py -- do the conclusions survive a leak in either accumulator?

The leak-robustness sweep (run_leak.sh) samples the leaks of the motion-evidence
accumulator (gamma_x) and of the response-signal accumulator (gamma_q) as two more
dimensions of a separate Latin-hypercube cube, uniform over [0, 1] /s.  Its draws are
screened with admissibility.py and classified exactly as the main sweep is (classify.py).

For every effect and model, the verdict is:
  reproduced     at least one admissible draw reproduces the observed magnitude
  unreproduced   every admissible draw is disjoint from the observed interval
  undecided      neither (no draw reproduces, some are uncertain)

The verdicts of the leak sweep are set against those of the main sweep (no leak).  The
two sweeps do not share their draws, so a difference in a count is not attributable to
the leaks alone; how each effect depends on the leaks is given within the leak sweep, by
the Spearman correlation, over the admissible draws, between the draw's effect and each
leak, and by the admissible counts in each third of the range of each leak.

Writes results/leak_summary.json and prints the comparison.
Usage:  python3 leak_summary.py [DUMP.jsonl]
"""
import json, os, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from effects import KEYS, NAMES, observed_per_subject, draw_effects
from classify import classify
from shardio import load as load_shards
from admissibility import reapply

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
LEAKS = os.path.join(HERE, 'results/leak/shard_*.jsonl')
OUT = os.path.join(HERE, 'results/leak_summary.json')
N_DRAWS = 1200
MODELS = ('model1', 'model2')
THIRDS = [(0.0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.0 + 1e-9)]


def verdict(r):
    if r['reproduces'] > 0:
        return 'reproduced'
    return 'unreproduced' if r['uncertain'] == 0 else 'undecided'


def summary(c):
    return dict(
        n_admissible={m: c['models'][m]['n_admissible'] for m in MODELS},
        verdicts={m: {k: verdict(c['models'][m][k]) for k in KEYS} for m in MODELS},
        counts={m: {k: {x: c['models'][m][k][x] for x in
                        ('n', 'min', 'max', 'reproduces', 'does_not', 'uncertain',
                         'sign_observed', 'sign_opposite', 'sign_uncertain')}
                    for k in KEYS} for m in MODELS})


def dependence(recs):
    """Per model: admissible counts by third of each leak, and, per effect, the Spearman
    correlation between the effect and each leak over the admissible draws."""
    out = {}
    for m in MODELS:
        adm = [r for r in recs if r[m]['admissible']]
        g = {name: np.array([r[name] for r in adm], float) for name in ('gamma_x', 'gamma_q')}
        allg = {name: np.array([r[name] for r in recs], float) for name in ('gamma_x', 'gamma_q')}
        thirds = {name: [dict(range=[lo, min(hi, 1.0)],
                              n_draws=int(((allg[name] >= lo) & (allg[name] < hi)).sum()),
                              n_admissible=int(((g[name] >= lo) & (g[name] < hi)).sum()))
                         for lo, hi in THIRDS] for name in g}
        est = [draw_effects(r[m]) for r in adm]
        corr = {}
        for k in KEYS:
            v = np.array([e[k] for e in est], float)
            ok = np.isfinite(v)
            corr[k] = {}
            for name in g:
                if ok.sum() >= 3:
                    rho, p = stats.spearmanr(g[name][ok], v[ok])
                    corr[k][name] = dict(rho=float(rho), p=float(p), n=int(ok.sum()))
                else:
                    corr[k][name] = dict(rho=None, p=None, n=int(ok.sum()))
        out[m] = dict(thirds=thirds, spearman=corr)
    return out


per_sub = observed_per_subject()
main_recs = [json.loads(l) for l in open(os.path.join(HERE, DUMP))]
leak_recs = [reapply(r) for r in load_shards(LEAKS)]
if len({r['draw'] for r in leak_recs}) < N_DRAWS:
    raise SystemExit(f"leak sweep: {len(leak_recs)} of {N_DRAWS} draws (run run_leak.sh)")

c_main, c_leak = classify(main_recs, per_sub), classify(leak_recs, per_sub)
base, leak = summary(c_main), summary(c_leak)
changes = [dict(model=m, effect=k, no_leak=base['verdicts'][m][k],
                leak_sampled=leak['verdicts'][m][k])
           for m in MODELS for k in KEYS
           if base['verdicts'][m][k] != leak['verdicts'][m][k]]
dep = dependence(leak_recs)
out = dict(source=DUMP, leak_source=os.path.relpath(LEAKS, HERE), n_draws=N_DRAWS,
           leak_range=[0.0, 1.0], observed=c_main['observed'], no_leak=base,
           leak_sampled=leak, changes=changes, dependence=dep)
json.dump(out, open(OUT, 'w'), indent=1)

print(f"{N_DRAWS} draws in each sweep. Admissible (Model 1 / Model 2): no leak "
      f"{base['n_admissible']['model1']}/{base['n_admissible']['model2']}, leaks sampled "
      f"{leak['n_admissible']['model1']}/{leak['n_admissible']['model2']}")
for m in MODELS:
    for name in ('gamma_x', 'gamma_q'):
        t = dep[m]['thirds'][name]
        print(f"  {m}, admissible by third of {name}: "
              + ', '.join(f"{x['n_admissible']}/{x['n_draws']}" for x in t))
lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
w = max(len(v) for v in lab.values())
print(f"\n{'effect':<{w}}  {'reproducing draws, no leak / leaks sampled':>44}"
      f"  {'Spearman rho with gamma_x, gamma_q':>36}")
for k in KEYS:
    cells, rhos = [], []
    for m in MODELS:
        a, b = base['counts'][m][k], leak['counts'][m][k]
        cells.append(f"M{m[-1]} {a['reproduces']}/{a['n']} | {b['reproduces']}/{b['n']}")
        s = dep[m]['spearman'][k]
        rhos.append(f"M{m[-1]} " + ' '.join(
            '  n/a' if s[x]['rho'] is None else f"{s[x]['rho']:+.2f}" for x in ('gamma_x', 'gamma_q')))
    print(f"{lab[k]:<{w}}  {'   '.join(cells):>44}  {'   '.join(rhos):>36}")
if changes:
    print("\nVerdicts that differ between the two sweeps:")
    for c in changes:
        print(f"  {c['model']}, {lab[c['effect']]}: {c['no_leak']} -> {c['leak_sampled']}")
else:
    print("\nNo verdict differs between the two sweeps, for either model.")
