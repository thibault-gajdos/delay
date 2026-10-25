"""
representative.py -- the parameter vector at which the single simulated trial of the
manuscript's Figure 2B is drawn.

Among the parameter vectors admissible in both models, the one closest to the
coordinate-wise median of the sampled parameters of that set, each coordinate standardized
by its interquartile range and taken in logarithm where it was sampled in logarithm.
beta is set to the median of the values admissible under Model 2, not to the value the
selected vector carries: requiring admissibility in both models selects smaller values of
beta, so the vector's own beta would not be typical of Model 2.  The rule involves the
parameters only, never a quantity the figure is used to show.

Writes results/representative.json.  Usage: python3 representative.py [DUMP.jsonl]
"""
import json, sys, numpy as np

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
OUT = 'results/representative.json'

recs = [json.loads(l) for l in open(DUMP)]

both = [r for r in recs if r['model1']['admissible'] and r['model2']['admissible']]
if not both:
    raise SystemExit('no draw is admissible in both models')

LOGDIMS = {'kappa', 'SDe', 'r', 'k', 'lam', 'sigq'}
SAMPLED = ['g', 'theta', 'kappa', 'SDe', 'muRS', 'r', 'k', 'lam', 'sigq']


def coords(r):
    v = []
    for name in SAMPLED:
        x = r['params'][name]
        v.append(np.log(x) if name in LOGDIMS else x)
    v.append(np.log(r['beta']))
    return np.array(v)


X = np.array([coords(r) for r in both])
mid = np.median(X, axis=0)
iqr = np.percentile(X, 75, axis=0) - np.percentile(X, 25, axis=0)
iqr[iqr <= 0] = 1.0
best = both[int(np.argmin((((X - mid) / iqr) ** 2).sum(axis=1)))]
bd = float((((coords(best) - mid) / iqr) ** 2).sum())

beta_m2 = float(np.median([r['beta'] for r in recs if r['model2']['admissible']]))
Delta = best['params']['delay'] + best['params']['LRS'] - best['params']['Lmotion']
alpha_m2 = beta_m2 / Delta

out = dict(source=DUMP, draw=best['draw'], distance=float(bd),
           params=best['params'], alpha=alpha_m2, beta=beta_m2,
           beta_of_draw=best['beta'],
           n_admissible_both=len(both),
           note='alpha/beta are the Model 2 admissible median, not the draw\'s own; '
                'see the docstring. beta_of_draw is what the draw carried.')
json.dump(out, open(OUT, 'w'), indent=1)
print(json.dumps(out, indent=1))
