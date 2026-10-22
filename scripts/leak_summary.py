"""leak_summary.py -- do the conclusions survive a leak in either accumulator?

The leaks of the motion-evidence accumulator (gamma_x) and of the response-signal
accumulator (gamma_q) are varied independently over gamma_x, gamma_q in {0, 0.25, 0.5, 1}
/s.  At every grid point the draws admissible in either model at zero leak are re-run
(chunk_leak.sh / run_leak.sh) and classified exactly as the main sweep is (classify.py).
Draws inadmissible at zero leak are not re-run, so a leak can remove admissible draws
but cannot add any.

For every effect and model, the verdict is:
  reproduced     at least one admissible draw reproduces the observed magnitude
  unreproduced   every admissible draw is disjoint from the observed interval
  undecided      neither (no draw reproduces, some are uncertain)

The zero-leak verdicts are computed on the same draws.  Writes results/leak_summary.json
and prints every verdict that differs from zero leak.

Usage:  python3 leak_summary.py [DUMP.jsonl]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from effect_panel import KEYS, NAMES, observed_per_subject
from classify import classify

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
G = ['0', '0.25', '0.5', '1']
OUT = os.path.join(HERE, 'results/leak_summary.json')


def load(path):
    return [json.loads(l) for l in open(path)]


def verdict(r):
    if r['reproduces'] > 0:
        return 'reproduced'
    return 'unreproduced' if r['uncertain'] == 0 else 'undecided'


idx = set(json.load(open(os.path.join(HERE, 'results/leak/adm_idx.json'))))
per_sub = observed_per_subject()
grid = {}
for gx in G:
    for gq in G:
        if gx == '0' and gq == '0':
            recs = [r for r in load(os.path.join(HERE, DUMP)) if r['draw'] in idx]
        else:
            d = os.path.join(HERE, f'results/leak/gx{gx}_gq{gq}')
            recs = [r for f in sorted(os.listdir(d)) if f.startswith('shard_')
                    for r in load(os.path.join(d, f))]
        if len(recs) < len(idx):
            raise SystemExit(f"gamma_x {gx}, gamma_q {gq}: {len(recs)} of {len(idx)} draws")
        c = classify(recs, per_sub)
        grid[(gx, gq)] = dict(
            n_admissible={m: c['models'][m]['n_admissible'] for m in ('model1', 'model2')},
            verdicts={m: {k: verdict(c['models'][m][k]) for k in KEYS}
                      for m in ('model1', 'model2')},
            counts={m: {k: {x: c['models'][m][k][x] for x in
                            ('reproduces', 'does_not', 'uncertain', 'sign_observed',
                             'sign_opposite', 'sign_uncertain')} for k in KEYS}
                    for m in ('model1', 'model2')})

base = grid[('0', '0')]
changes = []
for (gx, gq), g in grid.items():
    for m in ('model1', 'model2'):
        for k in KEYS:
            if g['verdicts'][m][k] != base['verdicts'][m][k]:
                changes.append(dict(gamma_x=float(gx), gamma_q=float(gq), model=m,
                                    effect=k, at_zero=base['verdicts'][m][k],
                                    at_leak=g['verdicts'][m][k]))
out = dict(source=DUMP, n_draws=len(idx), observed=classify([], per_sub)['observed'],
           grid=[dict(gamma_x=float(gx), gamma_q=float(gq), **g)
                 for (gx, gq), g in grid.items()],
           changes=changes)
json.dump(out, open(OUT, 'w'), indent=1)

print(f"{len(idx)} draws re-run at each leak. Admissible (Model 1 / Model 2):")
print('gamma_x \\ gamma_q ' + ''.join(f'{gq:>12}' for gq in G))
for gx in G:
    print(f'{gx:>17} ' + ''.join(
        f"{grid[(gx, gq)]['n_admissible']['model1']:>6}/{grid[(gx, gq)]['n_admissible']['model2']:<5}"
        for gq in G))
if changes:
    print("\nVerdicts that differ from zero leak:")
    for c in changes:
        print(f"  gamma_x {c['gamma_x']}, gamma_q {c['gamma_q']}, {c['model']}, "
              f"{NAMES[c['effect']].replace(chr(10), ' ')}: {c['at_zero']} -> {c['at_leak']}")
else:
    print("\nNo verdict differs from zero leak, for either model, at any grid point.")
