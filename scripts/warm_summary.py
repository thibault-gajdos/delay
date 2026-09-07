"""Summarise the stationary-start sweep against the observed delay-period statistics.

Run after run_warm_sweep.sh:  python3 warm_summary.py
Reads results/warm/shard_*.jsonl, writes results/sweep_screened_warm.jsonl and
results/warm_summary.json.

Observed targets are the raw-rate contrasts on the 20 delay bins of
PMA_rate_delay_JASP.csv, in the same units the sweep reports.  Their F-statistics
reproduce Table 4 of the manuscript exactly (35.238 and 18.891), which is what
certifies that the contrast definition matches the reported analysis.
"""
import json, glob, numpy as np

OBS = dict(rate=0.7091, trend=0.5604, curv=0.2099)      # s^-1, s^-2, s^-3

rows = []
for f in sorted(glob.glob("results/warm/shard_*.jsonl")):
    rows += [json.loads(l) for l in open(f)]
rows.sort(key=lambda r: r["draw"])
assert len({r["draw"] for r in rows}) == len(rows), "duplicate draw indices"
with open("results/sweep_screened_warm.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")

def cell(a, m, k):
    return np.asarray(a[m]["raw"][k], float)

out = {"n_draws": len(rows), "observed": OBS, "models": {}}
for m in ("model1", "model2"):
    adm = [r for r in rows if r[m]["admissible"]]
    d = {"admissible": len(adm)}
    for k in ("rate", "trend", "curv"):
        v = np.array([np.nanmean(cell(a, m, k)) for a in adm])
        v = v[np.isfinite(v)]
        d[k] = dict(median=float(np.median(v)), p5=float(np.percentile(v, 5)),
                    p95=float(np.percentile(v, 95)), max=float(v.max()),
                    frac_ge_obs=float((v >= OBS[k]).mean()))
    # the ramp at a rate within a factor two of the observed one
    rate = np.array([np.nanmean(cell(a, m, "rate")) for a in adm])
    tr = np.array([np.nanmean(cell(a, m, "trend")) for a in adm])
    near = np.isfinite(rate) & np.isfinite(tr) & (rate > OBS["rate"]/2) & (rate < OBS["rate"]*2)
    d["trend_at_observed_rate"] = dict(n=int(near.sum()),
                                       max=float(tr[near].max()) if near.any() else None)
    # trend significantly positive in every condition (the R1_v3 diagnostic)
    k = sum(1 for a in adm
            if np.all(np.isfinite(cell(a, m, "trend"))) and np.all(np.isfinite(cell(a, m, "trend_se")))
            and np.all(cell(a, m, "trend") - 1.96 * cell(a, m, "trend_se") > 0))
    d["all_six_positive"] = dict(n=k, frac=k/len(adm))
    # level effects
    R = np.array([cell(a, m, "rate") for a in adm])           # cells: 0.3s x3, then 2s x3
    for nm, v in (("L1_duration", R[:, 3:].mean(1) - R[:, :3].mean(1)),
                  ("L2_coherence", R[:, [2, 5]].mean(1) - R[:, [0, 3]].mean(1)),
                  ("L3_interaction", (R[:, 5] - R[:, 3]) - (R[:, 2] - R[:, 0]))):
        v = v[np.isfinite(v)]
        d[nm] = dict(median=float(np.median(v)), frac_positive=float((v > 0).mean()))
    out["models"][m] = d

json.dump(out, open("results/warm_summary.json", "w"), indent=2)
print(json.dumps(out, indent=2))
