"""g2_uncertainty.py -- how big is the Model 1 / Model 2 G^2 difference compared with
the Monte-Carlo noise in G^2 itself, and where in the objective does it come from?

The calibration uses common random numbers (fixed seeds 2000+k), so the objective is a
deterministic function of the parameters and the optimizer sees a smooth surface.  But
the VALUE of G^2 at the optimum is still a draw from a distribution over seeds, and the
difference between the two models' values inherits that noise.  This script re-evaluates
both models at their own fitted parameters over many seed offsets, and decomposes the
difference by response family and by cell.

Appends one JSON record per replicate to results/g2_uncertainty.jsonl (resumable).
Usage: python3 g2_uncertainty.py N_REPLICATES [N_TRIALS]
"""
import numpy as np, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_behav_q as C

N_ADD = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 25000     # matches the final stage
OUT = 'results/g2_uncertainty.jsonl'

P1 = json.load(open('results/calibq_m1.json'))
P2 = json.load(open('results/calibq_m2.json'))


def objective_seeded(p, model, n, seed_base):
    """C.objective with a shiftable seed base, so the same parameters can be re-evaluated
    under independent noise.  Mirrors C.objective exactly otherwise."""
    total, parts = 0.0, {}
    for k, cell in enumerate(C.TGT['cells']):
        out = C.simulate_behav(p, cell['coh'], cell['dur'], n,
                               model=model, seed=seed_base + k)
        rt_c, rt_i, n_prem = C.behav_rt_split(out)
        g2, g2_parts = C.cell_g2(cell, rt_c, rt_i, n_prem)
        parts[cell['label']] = g2_parts
        total += g2
    return total, parts


done = sum(1 for _ in open(OUT)) if os.path.exists(OUT) else 0
with open(OUT, 'a') as fh:
    for r in range(done, done + N_ADD):
        t0 = time.time()
        base = 2000 + 7919 * (r + 1)        # disjoint from the calibration's own seeds
        g1, p1 = objective_seeded(P1, 1, NTR, base)
        g2, p2 = objective_seeded(P2, 2, NTR, base)
        rec = dict(rep=r, seed_base=base, n=NTR, g2_m1=g1, g2_m2=g2, delta=g1 - g2,
                   parts_m1=p1, parts_m2=p2)
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        print(f"rep {r}: G2(M1)={g1:.3f}  G2(M2)={g2:.3f}  delta={g1-g2:+.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)
