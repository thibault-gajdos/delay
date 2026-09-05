"""leak_summary.py -- reduce the leak reruns to the table quoted in the response letter.

Reads results/sweep_screened.jsonl (gamma = 0) and results/leak/leak_g{025,050,100}.jsonl,
and writes results/leak_summary.json.  The raw leak dumps are regenerable and untracked;
this summary is what the letter cites.
"""
import json, os, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def load(path):
    return {json.loads(l)['draw']: json.loads(l) for l in open(os.path.join(HERE, path))}

base = load('results/sweep_screened.jsonl')
runs = {0.0: base}
for tag, g in (('025', 0.25), ('050', 0.50), ('100', 1.00)):
    runs[g] = load(f'results/leak/leak_g{tag}.jsonl')

adm0 = {m: {d for d, r in base.items() if r[m]['admissible']} for m in ('model1', 'model2')}
out = {'n_admissible_gamma0': {m: len(v) for m, v in adm0.items()}, 'rows': []}
for g in sorted(runs):
    R = runs[g]
    row = dict(gamma=g)
    row['m1_admissible'] = sum(1 for d in adm0['model1'] if d in R and R[d]['model1']['admissible'])
    row['m2_admissible'] = sum(1 for d in adm0['model2'] if d in R and R[d]['model2']['admissible'])
    tot = sum(1 for d in adm0['model2'] if d in R)
    row['n_evaluated'] = tot
    row['rt_coh_supported'] = sum(R[d]['model2']['rt']['coh_supported'] for d in adm0['model2'] if d in R)
    row['rt_dur_supported'] = sum(R[d]['model2']['rt']['dur_supported'] for d in adm0['model2'] if d in R)
    row['m2_all_six_slopes_positive'] = sum(1 for d in adm0['model2'] if d in R and R[d]['model2']['slope']['all_positive'])
    row['m1_all_six_slopes_positive'] = sum(1 for d in adm0['model1'] if d in R and R[d]['model1']['slope']['all_positive'])
    out['rows'].append(row)

R = runs[1.0]
lost = [d for d in adm0['model2'] if d in R and not R[d]['model2']['admissible']]
kept = [d for d in adm0['model2'] if d in R and R[d]['model2']['admissible']]
ratio = lambda ds: float(np.median([base[d]['params']['k'] * 0.5 / base[d]['params']['theta'] for d in ds]))
out['dropout_at_gamma_1'] = dict(n=len(lost), median_kc_over_theta_lost=ratio(lost),
                                 median_kc_over_theta_kept=ratio(kept),
                                 min_kc_over_theta_over_admissible=float(min(
                                     base[d]['params']['k'] * 0.5 / base[d]['params']['theta']
                                     for d in adm0['model1'] | adm0['model2'])))
with open(os.path.join(HERE, 'results/leak_summary.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
