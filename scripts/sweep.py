"""
sweep.py -- robustness of the directional predictions around the calibrated point.

Reviewer 2 asked that the qualitative predictions be shown to hold across a range of
parameters rather than at one hand-picked setting.  Each free parameter is multiplied
by 0.5 and by 2 in turn, everything else held at the calibrated value, and every
directional effect is re-evaluated.  A leak common to both accumulators is added as a
further variant (Reviewer 2's second point), at the rates the manuscript already uses.

The effects tested are the ones catalogued in the provenance document, restricted to
those the simulations bear on.  For each, the table records whether the sign matches
the data.
"""
import numpy as np, json, sys, itertools, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import PARAM_NAMES
from calibrate import model_stats, NTR
from delaymodel.targets import LAB

EFFECTS = [
    ('R1  PMET, duration',        lambda M: M['pmet'][:3].mean() - M['pmet'][3:].mean(), +1),
    ('R2  PMET, coherence',       lambda M: M['pmet'][[0, 3]].mean() - M['pmet'][[2, 5]].mean(), +1),
    ('R4  MET, no duration effect', lambda M: -abs(M['met'][:3].mean() - M['met'][3:].mean()) + 0.005, +1),
    ('R12 delay rate, duration',  lambda M: M['dly'][3:].mean() - M['dly'][:3].mean(), +1),
    ('R13 delay rate, coherence', lambda M: M['dly'][[2, 5]].mean() - M['dly'][[0, 3]].mean(), +1),
    ('R14 delay rate, time',      lambda M: M['slope'][np.isfinite(M['slope'])].mean(), +1),
    ('R15 delay rate, dur x coh', lambda M: (M['dly'][5] - M['dly'][3]) - (M['dly'][2] - M['dly'][0]), +1),
    ('R16 delay rate, dur x time', lambda M: M['slope'][3:].mean() - M['slope'][:3].mean(), -1),
    ('acc, duration at 12%',      lambda M: M['acc'][4] - M['acc'][1], +1),
    ('premature, duration',       lambda M: M['prem'][3:].mean() - M['prem'][:3].mean(), +1),
]


def evaluate(p, model=2, gamma=0.0, n=NTR):
    M = model_stats(p, n=n, model=model, gamma=gamma)
    return {name: f(M) for name, f, _ in EFFECTS}, M


def run(p0, model=2, n=NTR, factors=(0.5, 2.0)):
    rows = {}
    base, _ = evaluate(p0, model, n=n)
    rows['calibrated'] = base
    for k in PARAM_NAMES:
        for f in factors:
            q = dict(p0); q[k] = q[k] * f
            rows[f'{k} x{f}'] = evaluate(q, model, n=n)[0]
    for g in (0.25, 0.50, 1.00):
        rows[f'leak {g}'] = evaluate(p0, model, gamma=g, n=n)[0]
    return rows


def table(rows):
    names = [e[0] for e in EFFECTS]; sgn = {e[0]: e[2] for e in EFFECTS}
    hdr = f"{'setting':>16} " + " ".join(f"{n.split()[0]:>5}" for n in names)
    print(hdr); print('-' * len(hdr))
    ok = {n: 0 for n in names}; tot = 0
    for lab, r in rows.items():
        tot += 1
        cells = []
        for n in names:
            good = np.sign(r[n]) == sgn[n]
            ok[n] += bool(good)
            cells.append(' ok  ' if good else ' --  ')
        print(f"{lab:>16} " + "".join(cells))
    print('-' * len(hdr))
    print(f"{'holds in':>16} " + " ".join(f"{ok[n]:>2}/{tot:<2}" for n in names))


if __name__ == '__main__':
    p0 = json.load(open('results/fit_model2.json'))
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rows = run(p0, model=2, n=n)
    print("\nMODEL 2 (with temporal preparation)")
    table(rows)
    json.dump({k: {kk: float(vv) for kk, vv in v.items()} for k, v in rows.items()},
              open('results/sweep_model2.json', 'w'), indent=1)
    try:
        p1 = json.load(open('results/fit_model1.json'))
        rows1 = run(p1, model=1, n=n)
        print("\nMODEL 1 (no temporal preparation, recalibrated in its own right)")
        table(rows1)
        json.dump({k: {kk: float(vv) for kk, vv in v.items()} for k, v in rows1.items()},
                  open('results/sweep_model1.json', 'w'), indent=1)
    except FileNotFoundError:
        print("\n(Model 1 fit not found; run calibrate_m1.py first)")
