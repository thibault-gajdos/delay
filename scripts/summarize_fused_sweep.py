"""
summarize_fused_sweep.py -- distill a fused_sweep.py .jsonl dump into a small,
version-controlled provenance summary (the raw dump itself is regenerable and
gitignored, per scripts/README.md's "results the manuscript cites, not scratch"
convention -- this is the citable number, not the scratch behind it).

Usage: python3 summarize_fused_sweep.py IN.jsonl OUT.json
"""
import json, sys, numpy as np

IN = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None

recs = [json.loads(l) for l in open(IN)]
n = len(recs)


def discriminates(m):
    if not m['admissible']:
        return False
    return bool(m['rt']['coh_supported'] and m['rt']['dur_supported']
                and m['slope']['all_positive'])


summary = dict(source=IN, n_draws=n)
for tag, key in [('model1', 'model1'), ('model2', 'model2')]:
    ms = [r[key] for r in recs]
    adm = [m for m in ms if m['admissible']]
    disc = [m for m in ms if discriminates(m)]
    n_adm = len(adm)
    row = dict(n_admissible=n_adm,
               n_discriminates=len(disc),
               n_coh_supported=sum(1 for m in adm if m['rt']['coh_supported']),
               n_dur_supported=sum(1 for m in adm if m['rt']['dur_supported']),
               n_slope_all_positive=sum(1 for m in adm if m['slope']['all_positive']))
    if disc:
        minsl = np.array([m['slope']['min_slope'] for m in disc])
        maxsl = np.array([m['slope']['max_slope'] for m in disc])
        row['weakest_cell_slope'] = dict(min=float(minsl.min()), median=float(np.median(minsl)),
                                          max=float(minsl.max()))
        row['strongest_cell_slope'] = dict(min=float(maxsl.min()), median=float(np.median(maxsl)),
                                            max=float(maxsl.max()))
    met = np.array([m['met']['range_ms'] for m in adm])
    row['met_range_ms'] = dict(median=float(np.median(met)), p95=float(np.percentile(met, 95)),
                                max=float(met.max()))
    summary[tag] = row

text = json.dumps(summary, indent=1)
print(text)
if OUT:
    with open(OUT, 'w') as f:
        f.write(text)
