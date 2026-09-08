"""leak_summary.py -- reduce the leak reruns to the table quoted in the response letter.

Reads the zero-leak sweep dump and the three leak dumps beside it, and writes
results/leak_summary.json.  The raw leak dumps are regenerable and untracked; this
summary is what the letter cites.

  python3 leak_summary.py [DUMP.jsonl]

Both criteria are reported per model -- the raw-rate TREND, which is what the Results
cite, and the log-rate SLOPE, kept because it is what the earlier runs recorded and the
two should be seen to agree.
"""
import json, os, sys, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
LEAKDIR = 'results/leak'
OUT = 'results/leak_summary.json'


def load(path):
    return {json.loads(l)['draw']: json.loads(l) for l in open(os.path.join(HERE, path))}


base = load(DUMP)
runs = {0.0: base}
for tag, g in (('025', 0.25), ('050', 0.50), ('100', 1.00)):
    runs[g] = load(f'{LEAKDIR}/leak_g{tag}.jsonl')

adm0 = {m: {d for d, r in base.items() if r[m]['admissible']} for m in ('model1', 'model2')}
out = {'source': DUMP, 'leak_dir': LEAKDIR,
       'n_admissible_gamma0': {m: len(v) for m, v in adm0.items()}, 'rows': []}


def crit(rec, model, key):
    """`all_positive` under the named criterion, or None where the record predates it.

    The zero-leak row comes from the sweep dump, whose records hold the raw per-condition
    values rather than a criterion block, so it is recomputed here from those values --
    the same 95% lower bound leak_sweep applies."""
    blk = rec[model].get(key)
    if isinstance(blk, dict):
        return bool(blk['all_positive'])
    raw = rec[model].get('raw')
    if raw is None or key not in raw:
        return None
    v = np.asarray(raw[key], float)
    e = np.asarray(raw[f'{key}_se'], float)
    return bool(np.isfinite(v - e).all() and (v - 1.96 * e > 0).all())


for g in sorted(runs):
    R = runs[g]
    row = dict(gamma=g)
    for m, tag in (('model1', 'm1'), ('model2', 'm2')):
        row[f'{tag}_admissible'] = sum(1 for d in adm0[m] if d in R and R[d][m]['admissible'])
        for key in ('trend', 'slope'):
            vals = [crit(R[d], m, key) for d in adm0[m] if d in R]
            row[f'{tag}_all_six_{key}s_positive'] = (
                None if any(v is None for v in vals) else sum(vals))
    row['n_evaluated'] = sum(1 for d in adm0['model2'] if d in R)
    row['rt_coh_supported'] = sum(R[d]['model2']['rt']['coh_supported']
                                  for d in adm0['model2'] if d in R)
    row['rt_dur_supported'] = sum(R[d]['model2']['rt']['dur_supported']
                                  for d in adm0['model2'] if d in R)
    out['rows'].append(row)

R = runs[1.0]
lost = [d for d in adm0['model2'] if d in R and not R[d]['model2']['admissible']]
kept = [d for d in adm0['model2'] if d in R and R[d]['model2']['admissible']]
ratio = lambda ds: float(np.median([base[d]['params']['k'] * 0.5 / base[d]['params']['theta']
                                    for d in ds])) if ds else float('nan')
out['dropout_at_gamma_1'] = dict(
    n=len(lost), median_kc_over_theta_lost=ratio(lost), median_kc_over_theta_kept=ratio(kept),
    min_kc_over_theta_over_admissible=float(min(
        base[d]['params']['k'] * 0.5 / base[d]['params']['theta']
        for d in adm0['model1'] | adm0['model2'])))

with open(os.path.join(HERE, OUT), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
