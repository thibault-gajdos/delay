"""
representative.py -- picks the parameter vector at which the illustrative figures are
drawn, now that nothing is fitted.

The figures that show simulated dynamics (the manuscript's Figure 1 and Figure 2, and the
simulated rate profile of Figure S2) need one point in parameter space. It used to be a
fitted vector. With the fits removed it is instead a REPRESENTATIVE ADMISSIBLE DRAW of the
sweep, chosen by a rule that is stated in one line and involves no contact with the
quantities the figures are used to discuss.

Two properties are wanted, and one rules out the obvious choice.

  * It must be an actual draw, not a coordinate-wise median of the ranges. The median of
    the box is not admissible: with the ranges of Table 1 it sits at a standardized
    distance to the gate near 4, where exp(-a^2/2) leaves nothing to cross the gate and
    the model produces no muscle activity at all.

  * It must be the SAME vector in both models, so that the two panels of Figure 1 differ
    only in beta -- which is the whole point of the comparison. The sweep evaluates every
    draw in both models, so such draws exist; the draw is required to be admissible in
    both.

  * beta is set separately, to the median of the values admissible under Model 2 (0.0183),
    and not to the value the selected draw happens to carry. The reason is a selection
    effect worth recording: requiring admissibility in BOTH models biases beta downward
    (median 0.0075 there against 0.0183 among Model 2's own admissible draws), because a
    large ramp is what produces premature responses and the screen removes those. Taking
    beta from the draw would therefore illustrate Model 2 at a ramp smaller than is typical
    of Model 2, for a reason that has nothing to do with Model 2. Since beta is the ONLY
    parameter in which the two models differ, fixing it separately keeps the comparison
    exact rather than confounding it.

The rule: among draws admissible in both models, take the one closest to the
coordinate-wise median of the sampled parameters of that same set, each coordinate
standardized by its own interquartile range and taken in logarithm where it was sampled in
logarithm. In words, the most typical admissible parameter setting. The rule touches only
the parameters -- never the delay-period rate, the slope, or anything else the figures are
used to discuss -- so the illustration cannot be accused of having been selected for the
behaviour it illustrates. A first version of this script selected instead the draw closest
to the median RATE and SLOPE, and is a good example of why that is the wrong instinct: it
returned a draw at beta = 0.0018, the very bottom of the sampled range, at which Model 2
has no ramp to show.

Writes results/representative.json.
Usage: python3 representative.py [DUMP.jsonl]
"""
import json, sys, numpy as np

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
OUT = 'results/representative.json'

recs = [json.loads(l) for l in open(DUMP)]


def stats(m):
    r = np.nanmean(np.asarray(m['raw']['rate'], float))
    s = np.nanmean(np.asarray(m['raw']['slope'], float))
    return r, s


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
