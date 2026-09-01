"""
pool_noise.py -- pools the per-draw slope replicates that a sharded effect_sizes.py run
prints, and decomposes their variance once over all draws.

Why this exists as a script rather than a shell one-liner: each shard of effect_sizes.py
computes its own variance decomposition over only its own handful of draws, and a
between-draw variance estimated from seven draws is worthless -- one shard reported a
Monte-Carlo share above 100%, which is impossible. The decomposition has to be done once,
on the pooled per-draw values, which the shards write to their logs.

  total = between-draw + within-draw(Monte-Carlo)

with the within-draw part estimated as the mean over draws of the variance across seeds,
and the between-draw part as the variance of the per-draw means corrected for the fact that
each mean is itself noisy (subtracting within/n_seeds).

Reads results/logs/es_sh*.log, writes results/noise_pooled.json.
Usage: python3 pool_noise.py
"""
import glob, re, json
import numpy as np

rows = {'model1': [], 'model2': []}
for f in sorted(glob.glob('results/logs/es_sh*.log')):
    for line in open(f):
        m = re.match(r'^(model[12]) draw (\d+): (.+)$', line.strip())
        if m:
            rows[m.group(1)].append((int(m.group(2)), [float(x) for x in m.group(3).split()]))

out = {}
for tag, rs in rows.items():
    seen = {}
    for draw, vals in rs:
        seen.setdefault(draw, vals)
    A = np.array(list(seen.values()), float)
    A = A[np.isfinite(A).all(axis=1)]
    if A.size == 0:
        continue
    n, k = A.shape
    means = A.mean(axis=1)
    within = A.var(axis=1, ddof=1).mean()
    total = means.var(ddof=1) + within * (1 - 1.0 / k)
    out[tag] = dict(n_draws=int(n), n_seeds=int(k),
                    sd_total=float(np.sqrt(total)), sd_mc=float(np.sqrt(within)),
                    sd_between=float(np.sqrt(max(total - within, 0.0))),
                    mc_share=float(within / total),
                    se_per_draw=float(np.sqrt(within / k)))

json.dump(out, open('results/noise_pooled.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
