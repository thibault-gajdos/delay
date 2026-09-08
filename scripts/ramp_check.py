"""ramp_check.py -- does Model 1 produce a delay-period ramp ANYWHERE in the sampled space?

The question has an exact answer -- Appendix Corollary A.3: Model 1's delay-period rate is
constant, so its trend is zero at every parameter setting -- so this script's job is to
check that the SIMULATION agrees with it once the Monte-Carlo error is small enough to see
the difference. The statistic is the raw-rate TREND, the one the Results report.

An existence question is badly served by raising the trial count uniformly across a
Latin-hypercube sweep: the statistic that answers it is a MAXIMUM over draws, and a
maximum over noisy estimates is biased upward, so precision belongs on the candidates
rather than spread evenly over the space.  Two stages therefore:

  stage 1  fused_sweep.py at moderate NTR, to FIND candidates
  stage 2  this script: re-simulate the best of them at high NTR, several independent
           seeds each, giving a trend with a Monte-Carlo standard error attached

Selects Model 1's best draws two ways -- the best among those producing a realistic
activation rate (the honest conditional question), and the best overall regardless of
rate (in case a ramp hides at an implausibly low rate) -- and reports both.

Appends one record per (draw, seed) to results/ramp_check.jsonl; resumable.

Each record carries `dump`, the sweep it was selected from, and records from any other dump
are ignored on load. Without that the file silently accumulated records from two different
Latin-hypercube cubes, in which the same `draw` index means two different parameter
vectors, and is not detectable by eye.
Usage: python3 ramp_check.py [N_SEEDS] [NTR] [DUMP.jsonl]
"""
import numpy as np, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
DUMP = sys.argv[3] if len(sys.argv) > 3 else 'results/sweep_screened.jsonl'
# fused_sweep.py parses sys.argv at import time; hide this script's arguments from it.
_argv, sys.argv = sys.argv, sys.argv[:1]
import fused_sweep as F
sys.argv = _argv
F.NTR = NTR
OUT = 'results/ramp_check.jsonl'
# Observed delay-period partial-activation statistics, all from the 132 subject x
# condition cells (22 participants x 6 cells, none empty; see README "The EMG criterion").
# OBS_RATE is the MEAN rate over the 1-3 s analysis window -- the quantity the simulation
# reports, and not the end-of-delay peak of the last 0.1 s bin, which is a different
# quantity (1.446 /s) and would overstate both models' shortfall by about a factor of two.
OBS_RATE, OBS_RATE_MIN, OBS_RATE_MAX = 0.709, 0.085, 2.785
OBS_SLOPE_LO = 0.659

recs = [json.loads(l) for l in open(DUMP)]
cand = []
for r in recs:
    m = r['model1']
    if not m['admissible']:
        continue
    sl = np.asarray(m['raw']['trend'], float); rt = np.asarray(m['raw']['rate'], float)
    if not np.isfinite(sl).any():
        continue
    cand.append(dict(draw=r['draw'], params=r['params'],
                     trend=float(np.nanmean(sl)), rate=float(np.nanmean(rt))))
withrate = sorted([c for c in cand if c['rate'] >= 0.25 * OBS_RATE],
                  key=lambda c: -c['trend'])[:5]
overall = sorted(cand, key=lambda c: -c['trend'])[:5]
sel, seen = [], set()
for c, tag in [(c, 'rate>=25%') for c in withrate] + [(c, 'any rate') for c in overall]:
    if c['draw'] in seen:
        continue
    seen.add(c['draw']); c = dict(c); c['why'] = tag; sel.append(c)

done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        d = json.loads(l)
        if d.get('dump') != DUMP:
            continue
        done.add((d['draw'], d['seed']))

with open(OUT, 'a') as fh:
    for c in sel:
        for s in range(N_SEEDS):
            if (c['draw'], s) in done:
                continue
            t0 = time.time()
            res = F.run_model(c['params'], 0.0, 1, 777000 + 9000 * s + 13 * c['draw'])
            sl = np.asarray(res['trend'], float); rt = np.asarray(res['rate'], float)
            rec = dict(draw=c['draw'], seed=s, ntr=NTR, why=c['why'], dump=DUMP,
                       stage1_trend=c['trend'], stage1_rate=c['rate'],
                       trend=float(np.nanmean(sl)), rate=float(np.nanmean(rt)),
                       trend_cells=sl.tolist(), rate_cells=rt.tolist())
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            print(f"draw {c['draw']:3d} seed {s}: trend {rec['trend']:+.4f} "
                  f"rate {rec['rate']:.3f} ({time.time()-t0:.0f}s)", flush=True)

# The verdict, recomputed from the whole file each time so that an interrupted run still
# reports on what it has. The per-draw mean over seeds is the estimate; the spread over
# seeds is its Monte-Carlo error, and the question is only whether the estimate is
# distinguishable from the zero Corollary A.3 requires.
recs = [json.loads(l) for l in open(OUT) if json.loads(l).get('dump') == DUMP]
if recs:
    by = {}
    for r in recs:
        by.setdefault(r['draw'], []).append(r)
    print(f"\n{'draw':>5} {'rate':>7} {'seeds':>6} {'mean trend':>11} {'SE':>8}")
    worst = 0.0
    for d, rs in sorted(by.items(), key=lambda kv: -np.mean([x['trend'] for x in kv[1]])):
        t = np.array([x['trend'] for x in rs])
        se = t.std(ddof=1) / np.sqrt(len(t)) if len(t) > 1 else float('nan')
        worst = max(worst, float(t.mean()))
        print(f"{d:>5} {np.mean([x['rate'] for x in rs]):7.3f} {len(t):>6} "
              f"{t.mean():>+11.4f} {se:>8.4f}")
    print(f"\nlargest mean trend over {len(by)} candidates at NTR={NTR}: {worst:+.4f} /s^2, "
          f"against an observed +0.5604 and a predicted 0 (Appendix, Corollary A.3)")
