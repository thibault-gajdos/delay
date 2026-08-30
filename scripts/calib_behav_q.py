"""Behavioural calibration, Ratcliff quantile / G^2 version (replaces calib_behav.py).

For each of the 6 coherence x duration cells, three response-type families are scored:
premature, correct, incorrect. Correct and incorrect are further split by the 5 data
quantiles (.1/.3/.5/.7/.9, from targets_q.json), partitioning each response type into 6
bins holding, by construction, 10/20/20/20/20/10% of the DATA's trials of that type.
Premature is a single proportion with no RT substructure. The model is scored by how far
the SIMULATED RTs, passed through those same (data-fixed) bin edges, land from that
reference -- together with how far the simulated premature/correct/incorrect proportions
are from the data's -- via the G^2 statistic

    G^2 = 2 * sum_bins  n_obs * ln(n_obs / n_pred)

summed across all three families and all 6 cells (see g2_term/cell_g2 below -- every
family goes through the same g2_term, which is what makes the sum a valid, additive
deviance rather than an ad hoc combination).

n_obs and n_pred are PSEUDO-counts, not raw trial counts: the per-cell sample size is
fixed at n_subj = 22 (the two-stage unit -- one participant, one vote), not the pooled
trial count, so that a cell's weight in the objective does not depend on how many trials
happened to be collected in it. This keeps the two-stage-aggregation principle used
throughout the rest of the calibration; it means G^2 is used here purely as a scalar
discrepancy for calibration, not as a chi-square test with a meaningful p-value or df.

Premature responses were reintroduced into the objective 2026-08-30 (previously excluded
entirely, alongside <0.2s and no-response, going back to when this module replaced
calib_behav.py): theta, g, r and kappa -- the model's least independently-constrained
parameters -- are exactly the ones that set the premature rate, and until now nothing in
this objective carried any information about it at all, only a loose 10% ceiling three
steps downstream in the parameter-sweep's admissibility screen. <0.2s and
no-response-by-deadline trials remain excluded entirely, not scored, not counted in any
denominator -- the historical finding that no-response-by-deadline is not reproducible at
any behaviourally-plausible setting is untouched by this change (and, checked while
implementing this: no-response-by-deadline does not occur as a distinct category at all in
this archive's raw trial records, 0/21120 -- see build_targets_q.py's docstring).
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


def g2_term(obs, pred):
    """2 * sum(obs * log(obs/pred)), both already N_SUBJ-scaled pseudo-counts. Every
    family below goes through this one function -- it, not anything family-specific, is
    what makes G^2 additive across independent families (and would make an EMG-derived
    family additive too, if one is ever added; see README, "Deferred")."""
    obs = np.maximum(obs, FLOOR * N_SUBJ)
    pred = np.maximum(pred, FLOOR * N_SUBJ)
    return 2.0 * np.sum(obs * np.log(obs / pred))


def cell_g2(cell, rt_cor, rt_incor, n_prem):
    """Three response-type families per cell: premature (a single proportion, no RT
    structure -- there is no quantile target for it), correct and incorrect (six
    RT-quantile bins each, exactly as before). Premature responses were reintroduced
    2026-08-30 (see this module's docstring): theta/g/r/kappa -- the parameters hardest
    to bound independently -- are exactly the ones that set how often the model responds
    before the RS, and until now nothing in this objective was sensitive to that at all.
    p_correct/p_incorrect in the target file keep their original meaning (accuracy
    conditional on a valid response); they are rescaled by (1-p_premature) here, at the
    point of use, to get the three ABSOLUTE proportions the G^2 actually needs."""
    n_c, n_i = rt_cor.size, rt_incor.size
    n_tot = max(n_c + n_i + n_prem, 1)
    p_sim_prem = n_prem / n_tot
    p_sim_c = n_c / n_tot; p_sim_i = n_i / n_tot

    p_obs_prem = cell['p_premature']
    p_obs_c = (1 - p_obs_prem) * cell['p_correct']
    p_obs_i = (1 - p_obs_prem) * cell['p_incorrect']

    families = {
        'premature': (np.array([p_obs_prem]) * N_SUBJ,
                      np.array([p_sim_prem]) * N_SUBJ),
        'correct':   (p_obs_c * BIN_P * N_SUBJ,
                      p_sim_c * bin_fracs(rt_cor, cell['correct_quantiles']) * N_SUBJ),
        'incorrect': (p_obs_i * BIN_P * N_SUBJ,
                      p_sim_i * bin_fracs(rt_incor, cell['incorrect_quantiles']) * N_SUBJ),
    }
    # Gibbs' inequality guarantees the CELL TOTAL is >= 0 (obs and pred are complete
    # distributions over all 13 categories, both summing to N_SUBJ) -- individual family
    # subtotals below are NOT so guaranteed and can and do come out negative (a family is
    # an arbitrary subset of those 13 categories, not a complete distribution on its own).
    # That is expected, not a bug; only the sum across all three is the meaningful number.
    parts = {name: g2_term(obs, pred) for name, (obs, pred) in families.items()}
    return sum(parts.values()), parts


def objective(p, model, n=20000, detail=False):
    total = 0.0; parts = {}
    for k, cell in enumerate(TGT['cells']):
        c, D = cell['coh'], cell['dur']
        out = simulate_behav(p, c, D, n, model=model, seed=2000 + k)
        rt_cor, rt_incor, n_prem = behav_rt_split(out)
        g2, g2_parts = cell_g2(cell, rt_cor, rt_incor, n_prem)
        parts[cell['label']] = g2_parts
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
