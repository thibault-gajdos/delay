"""augment_windows.py -- add the stimulus-, response-signal- and response-window summaries
to the sweep dump, for the draws the effect panel uses.

The sweep records, per condition, a level and a linear trend for the delay window and a
level for the stimulus window.  The Results report effects in the other three windows too
-- the coherence- and duration-dependence of the stimulus-period ramp, the fall of the rate
after the response signal and its duration-dependence, the level before the response --
and none of them could be drawn against the model.  This script fills that in without a
full re-sweep: only the draws admissible in either model matter to the panel, and they are
already listed in results/leak/adm_idx.json.

Each draw is re-simulated with the SAME seeds as the sweep (fused_sweep.run_model is
called with the same seed0), so the new per-condition summaries belong to the same
realisation as the ones already in the dump.  Output is one JSON line per draw with only
the new fields; effect_panel.py merges it onto the dump by draw index.  Resumable by line
count.

Usage:  python3 augment_windows.py [OUT.jsonl]
        SHARD_N / SHARD_I in the environment split the draw list across cores.
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results/window_effects.jsonl')
sys.argv = ['fused_sweep.py', '1200', '4000', '/dev/null', '0.0']
import fused_sweep as FS

FIELDS = ('stim_trend', 'stim_trend_se', 'rs_rate', 'rs_rate_se',
          'rs_trend', 'rs_trend_se', 'resp_rate', 'resp_rate_se',
          'resp_trend', 'resp_trend_se')
IDX = json.load(open(os.path.join(HERE, 'results/leak/adm_idx.json')))
sh_n = int(os.environ.get('SHARD_N', '1'))
sh_i = int(os.environ.get('SHARD_I', '0'))
IDX = [d for j, d in enumerate(IDX) if j % sh_n == sh_i]

rng = np.random.default_rng(FS.SEED0)
cube = FS.lhs(FS.N_DRAWS, FS.DIMS, rng)
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
done = sum(1 for _ in open(OUT)) if os.path.exists(OUT) else 0
print(f"shard {sh_i}/{sh_n}: {len(IDX)} draws, resuming at {done}", flush=True)

with open(OUT, 'a') as f:
    for i in IDX[done:]:
        t0 = time.time()
        p, bf = FS.to_param(cube[i])
        beta2 = np.exp(FS.LO['log_beta'] + bf * (FS.HI['log_beta'] - FS.LO['log_beta']))
        Delta = (FS.FIXED['delay'] + FS.FIXED['LRS']) - FS.FIXED['Lmotion']
        r1 = FS.run_model(p, 0.0, 1, FS.SEED0 + 1000 * i)
        r2 = FS.run_model(p, beta2 / Delta, 2, FS.SEED0 + 1000 * i + 500)
        rec = dict(draw=i)
        for tag, r in (('model1', r1), ('model2', r2)):
            rec[tag] = {k: np.asarray(r[k], float).tolist() for k in FIELDS}
        f.write(json.dumps(rec) + "\n")
        f.flush()
        print(f"draw {i} ({time.time()-t0:.1f}s)", flush=True)
print("SHARD DONE", flush=True)
