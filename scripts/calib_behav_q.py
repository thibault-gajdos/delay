"""Behavioural calibration, Ratcliff quantile / G^2 version (replaces calib_behav.py).

For each of the 6 coherence x duration cells, and separately for correct and incorrect
responses, the 5 data quantiles (.1/.3/.5/.7/.9, from targets_q.json) partition that
response type into 6 bins holding, by construction, 10/20/20/20/20/10% of the DATA's
trials of that type. The model is scored by how far the SIMULATED RTs, passed through
those same (data-fixed) bin edges, land from that reference -- together with how far the
simulated accuracy is from the data's -- via the G^2 statistic

    G^2 = 2 * sum_bins  n_obs * ln(n_obs / n_pred)

n_obs and n_pred are PSEUDO-counts, not raw trial counts: the per-cell sample size is
fixed at n_subj = 22 (the two-stage unit -- one participant, one vote), not the pooled
trial count, so that a cell's weight in the objective does not depend on how many trials
happened to be collected in it. This keeps the two-stage-aggregation principle used
throughout the rest of the calibration; it means G^2 is used here purely as a scalar
discrepancy for calibration, not as a chi-square test with a meaningful p-value or df.

Premature, <0.2s and no-response trials are excluded entirely (not scored): p(correct) +
p(incorrect) = 1 over valid trials only, both in the data and in the simulation.
"""
import numpy as np, json, os, sys, time
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.behav_sim import simulate_behav, behav_rt_split

N_SUBJ = 22.0          # two-stage pseudo sample size per cell
FLOOR = 1e-4            # minimum bin/response-type probability, avoids log(0) / G^2 = inf
BIN_P = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1])

TGT = json.load(open('results/targets_q.json'))
CELLS = [(c['coh'], c['dur']) for c in TGT['cells']]
LAB = [c['label'] for c in TGT['cells']]


def bin_fracs(rt, edges):
    """Fraction of `rt` falling in each of the 6 bins defined by the 5 `edges`."""
    if rt.size == 0:
        return np.full(6, FLOOR)
    idx = np.searchsorted(edges, rt, side='right')
    cnt = np.bincount(idx, minlength=6)[:6].astype(float)
    return np.maximum(cnt / rt.size, FLOOR)


def cell_g2(cell, rt_cor, rt_incor):
    p_obs_c = cell['p_correct']; p_obs_i = cell['p_incorrect']
    n_c, n_i = rt_cor.size, rt_incor.size
    n_tot = max(n_c + n_i, 1)
    p_sim_c = max(n_c / n_tot, FLOOR); p_sim_i = max(n_i / n_tot, FLOOR)

    obs = np.concatenate([p_obs_c * BIN_P, p_obs_i * BIN_P]) * N_SUBJ
    pred_frac = np.concatenate([bin_fracs(rt_cor, cell['correct_quantiles']),
                                 bin_fracs(rt_incor, cell['incorrect_quantiles'])])
    pred = np.concatenate([np.full(6, p_sim_c), np.full(6, p_sim_i)]) * pred_frac * N_SUBJ
    pred = np.maximum(pred, FLOOR * N_SUBJ)
    obs = np.maximum(obs, FLOOR * N_SUBJ)
    return 2.0 * np.sum(obs * np.log(obs / pred))


def objective(p, model, n=20000, detail=False):
    total = 0.0; parts = {}
    for k, cell in enumerate(TGT['cells']):
        c, D = cell['coh'], cell['dur']
        out = simulate_behav(p, c, D, n, model=model, seed=2000 + k)
        rt_cor, rt_incor = behav_rt_split(out)
        g2 = cell_g2(cell, rt_cor, rt_incor)
        parts[cell['label']] = g2
        total += g2
    return (total, parts) if detail else total


STAGES = [dict(n=6000, scale=0.10, maxfev=90),
          dict(n=12000, scale=0.05, maxfev=70),
          dict(n=25000, scale=0.02, maxfev=50)]


def run(model, start, tag):
    names = list(start.keys())
    ckpt = f'results/ckptq_{tag}.json'
    stage0 = 0
    if os.path.exists(ckpt):
        st = json.load(open(ckpt))
        start = st['p']; stage0 = st.get('stage', 0)
    x0 = np.log([start[k] for k in names])
    t0 = time.time()
    p = dict(zip(names, np.exp(x0))); fun = None
    for stage in range(stage0, len(STAGES)):
        cfg = STAGES[stage]

        def f(v, n=cfg['n']):
            return objective(dict(zip(names, np.exp(v))), model, n=n)

        sim0 = np.vstack([x0] + [x0 + cfg['scale'] * np.eye(len(x0))[i]
                                  for i in range(len(x0))])
        res = minimize(f, x0, method='Nelder-Mead',
                        options=dict(maxfev=cfg['maxfev'], xatol=1e-4, fatol=1e-3,
                                     adaptive=True, initial_simplex=sim0))
        x0 = res.x; fun = res.fun
        p = dict(zip(names, np.exp(x0)))
        json.dump(dict(p=p, obj=float(res.fun), stage=stage + 1, nfev=int(res.nfev)),
                   open(ckpt, 'w'), indent=1)
        elapsed = (time.time() - t0) / 60
        with open(f'results/logs/logq_{tag}.txt', 'a') as fh:
            fh.write(f"stage {stage} n={cfg['n']} scale={cfg['scale']}: obj={res.fun:.3f} "
                     f"nfev={res.nfev} elapsed={elapsed:.1f}min\n")
    return p, fun


if __name__ == '__main__':
    tag = sys.argv[1]           # 'm1' or 'm2'
    model = 1 if tag == 'm1' else 2
    start = json.load(open(f'results/calibb_{tag}.json'))   # previous (MSE-based) calibration
    p, obj = run(model, start, tag)
    json.dump(p, open(f'results/calibq_{tag}.json', 'w'), indent=1)
    print(f"{tag}: G2={obj:.3f}")
