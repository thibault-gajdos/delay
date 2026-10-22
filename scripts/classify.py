"""
classify.py -- the model-data comparison of the Methods ("Model-data comparison").

For every effect of effect_panel.py and every admissible draw of each model, the draw's
estimate and its Monte-Carlo 95% interval (estimate +- 1.96 Monte-Carlo standard errors)
are set against the empirical 95% confidence interval of the group mean (Student t,
n - 1 degrees of freedom):

  reproduces       the Monte-Carlo interval lies entirely inside the empirical interval
  does not         the two intervals are disjoint
  uncertain        otherwise

A model reproduces an effect's magnitude if at least one admissible draw reproduces it; the
magnitude is unreproduced within the explored space only if every admissible draw does
not.  The sign of a draw's effect is the observed one, the opposite one, or uncertain,
according to whether its Monte-Carlo interval excludes zero and on which side.

Sufficiency of the trial count.  Uncertain draws matter only where they decide a verdict:
an effect that no draw reproduces, but that some draw is uncertain about, is neither
reproduced nor established as unreproduced.  Such effects are listed at the end; if there
are any, the whole sweep is to be rerun with more trials per condition.

Writes results/classification.json.  Usage: python3 classify.py [DUMP.jsonl]
"""
import json, sys, os
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effect_panel import (KEYS, GROUPS, NAMES, SRC, observed_per_subject, draw_effects,
                          effects_se)

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
OUT = 'results/classification.json'


def main():
    recs = [json.loads(l) for l in open(DUMP)]
    per_sub = observed_per_subject()
    out = dict(source=DUMP, n_draws=len(recs), observed={}, models={})
    for k in KEYS:
        x = np.array([p[k] for p in per_sub], float)
        x = x[np.isfinite(x)]
        n = x.size
        m, se = x.mean(), x.std(ddof=1) / np.sqrt(n)
        h = stats.t.ppf(0.975, n - 1) * se
        out['observed'][k] = dict(mean=float(m), lo=float(m - h), hi=float(m + h), n=n)
    undecided = []
    for tag in ('model1', 'model2'):
        adm = [r[tag] for r in recs if r[tag]['admissible']]
        est = [draw_effects(m) for m in adm]
        ses = [effects_se(m) for m in adm]
        res = dict(n_admissible=len(adm))
        for k in KEYS:
            v = np.array([e[k] for e in est], float)
            s = np.array([e[k] for e in ses], float)
            ok = np.isfinite(v) & np.isfinite(s)
            v, s = v[ok], s[ok]
            lo, hi = v - 1.96 * s, v + 1.96 * s
            o = out['observed'][k]
            inside = (lo >= o['lo']) & (hi <= o['hi'])
            disjoint = (hi < o['lo']) | (lo > o['hi'])
            unc = ~inside & ~disjoint
            sg = np.sign(o['mean'])
            same = (lo > 0) if sg > 0 else (hi < 0)
            opp = (hi < 0) if sg > 0 else (lo > 0)
            nt = int(v.size)
            res[k] = dict(
                n=nt, min=float(v.min()) if nt else None, max=float(v.max()) if nt else None,
                reproduces=int(inside.sum()), does_not=int(disjoint.sum()),
                uncertain=int(unc.sum()),
                sign_observed=int(same.sum()), sign_opposite=int(opp.sum()),
                sign_uncertain=int(nt - same.sum() - opp.sum()))
            if nt and inside.sum() == 0 and unc.sum() > 0:
                undecided.append((tag, k, int(unc.sum())))
        out['models'][tag] = res
    out['undecided'] = undecided
    json.dump(out, open(OUT, 'w'), indent=1)

    lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
    w = max(len(v) for v in lab.values())
    print(f"{out['n_draws']} draws; admissible: Model 1 {out['models']['model1']['n_admissible']},"
          f" Model 2 {out['models']['model2']['n_admissible']}")
    print(f"{'effect':{w}} {'observed 95% CI':>24}   {'M1 rep/not/unc':>15} {'M2 rep/not/unc':>15}")
    for win, ks in GROUPS:
        print(f"-- {win}")
        for k in ks:
            f = '%+.0f' if SRC[k][1] == 1000.0 else '%+.3f'
            o = out['observed'][k]
            ci = f"[{f % o['lo']}, {f % o['hi']}]"
            c = [f"{r['reproduces']}/{r['does_not']}/{r['uncertain']}"
                 for r in (out['models']['model1'][k], out['models']['model2'][k])]
            print(f"{lab[k]:{w}} {ci:>24}   {c[0]:>15} {c[1]:>15}")
    if undecided:
        print("\nUNDECIDED (no draw reproduces, some uncertain) -- more trials needed:")
        for tag, k, u in undecided:
            print(f"  {tag} {k}: {u} uncertain draws")
    else:
        print("\nNo verdict depends on an uncertain draw.")


if __name__ == '__main__':
    main()
