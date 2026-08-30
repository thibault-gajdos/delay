"""
calibrate.py -- calibration of the model to the electromyographic and behavioural data.

Procedure, in two stages.

STAGE A, closed form.  The identities of the note determine a starting point without
any search:
    omega   := sqrt(lambda kappa)   from the crossing rate,  nu = (omega/2pi) exp(-a^2/2)
    a                               from the count-only identity  abar = 1.43 sqrt(A)
    r/mu_RS := mean MET - L_motor
    (g-|m|-beta)/mu_RS := mean PMET - L_RS - 1/lambda
    sigma_Y := (g-|m|-beta)/a
    rho_r   := r/sigma_Y            checked against the premature-response rate through
                                    nu(g+r)/nu(g) = exp[-rho_r (a + rho_r/2)]
    alpha   := sigma_Y * slope / a  from the delay log-rate slope
Only ratios are determined, so sigma_x fixes the amplitude unit and one latency has to
be chosen; the search then relaxes everything.

STAGE B, search.  Nelder-Mead on all thirteen free parameters, minimising a normalised
squared error over nine blocks of six cells each: accuracy, mean and SD of premotor
execution time, mean and SD of motor execution time, delay-period partial-activation
rate, delay-period log-rate slope, stimulus-period rate, and the premature-response
rate.  Each block is divided by the between-participant standard error of the
corresponding target, so that a unit of error means one standard error.  Common random
numbers are used across evaluations so the objective is deterministic.
"""
import numpy as np, sys, time, json, os
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, rate_profile, PARAM_NAMES
from delaymodel.targets import CELLS, LAB, build, summary, poisson_slope

NTR = 200
TB = np.arange(20) * 0.1 + 0.05; TB = TB - TB.mean()


def model_stats(p, n=NTR, model=2, seed0=1000, gamma=0.0):
    out = {k: np.zeros(6) for k in ALL_BLOCKS}
    prof = np.zeros((6, 20))
    out['rt'] = np.zeros(6); out['rt_sd'] = np.zeros(6)
    for k, (c, D) in enumerate(CELLS):
        s = simulate(p, c, D, n, model, seed=seed0 + k, gamma=gamma)
        m = measures(s)
        v = m['valid']
        out['prem'][k] = m['premature'].mean()
        if v.sum() < 10:
            out['acc'][k] = np.nan; continue
        out['acc'][k] = m['acc'][v].mean()
        rv = m['rt'][v]
        out['rt'][k] = np.nanmean(rv)
        out['rt_sd'][k] = np.subtract(*np.nanpercentile(rv, [75, 25]))
        out['pmet'][k] = np.nanmean(m['pmet'][v]); out['met'][k] = np.nanmean(m['met'][v])
        pv = m['pmet'][v]; mv = m['met'][v]
        out['pmet_sd'][k] = np.subtract(*np.nanpercentile(pv, [75, 25]))
        out['met_sd'][k] = np.subtract(*np.nanpercentile(mv, [75, 25]))
        pr = rate_profile(m, D + 1, D + 3, 20); prof[k] = pr
        out['dly'][k] = pr.mean()
        out['stim'][k] = rate_profile(m, 0.0, 2.0, 20).mean()
        out['slope'][k] = poisson_slope(pr, TB)[1] if (pr > 0).sum() >= 4 else np.nan
    out['profile'] = prof
    return out


ALL_BLOCKS = ['acc', 'pmet', 'met', 'pmet_sd', 'met_sd', 'dly', 'stim', 'prem', 'slope']
BEHAV_BLOCKS = ['acc', 'rt', 'rt_sd', 'prem']      # what a behaviour-only study sees
BLOCKS = ALL_BLOCKS


def make_objective(E):
    S = summary(E)
    scale = {}
    for b in BLOCKS:
        key = {'dly': 'dly', 'stim': 'stim'}.get(b, b)
        arr = E[key]
        if arr.ndim == 3:
            arr = arr.mean(-1)
        se = np.nanstd(arr, 0, ddof=1) / np.sqrt(arr.shape[0])
        scale[b] = np.maximum(se, 1e-6)
    tgt = {b: S[b] for b in BLOCKS}

    def obj(p, n=NTR, model=2, detail=False, gamma=0.0):
        M = model_stats(p, n, model, gamma=gamma)
        tot = 0.0; parts = {}
        for b in BLOCKS:
            d = (M[b] - tgt[b]) / scale[b]
            d = d[np.isfinite(d)]
            parts[b] = float(np.mean(d ** 2)) if d.size else 100.0
            tot += parts[b]
        return (tot, parts, M) if detail else tot
    return obj, tgt, scale


def unpack(v):
    q = dict(zip(PARAM_NAMES, np.exp(v)))
    return q


def pack(p):
    return np.log([p[k] for k in PARAM_NAMES])


START = dict(k=1.10, lam=22.0, theta=0.085, g=0.150, kappa=11.8, SDe=0.038,
             alpha=0.0125, muRS=0.60, sigq=0.14, r=0.070, Lmotor=0.020,
             Lmotion=0.20, LRS=0.20)

LO = dict(k=0.3, lam=5.0, theta=0.01, g=0.02, kappa=1.0, SDe=0.005, alpha=0.001,
          muRS=0.05, sigq=0.01, r=0.005, Lmotor=0.005, Lmotion=0.05, LRS=0.05)
HI = dict(k=4.0, lam=60.0, theta=0.60, g=0.80, kappa=60.0, SDe=0.20, alpha=0.10,
          muRS=5.0, sigq=1.0, r=0.40, Lmotor=0.15, Lmotion=0.40, LRS=0.40)


def clip(p):
    return {k: float(min(max(p[k], LO[k]), HI[k])) for k in PARAM_NAMES}


if __name__ == '__main__':
    # chunked operation: each call runs a bounded number of evaluations, then
    # checkpoints, so the search can be resumed.
    import os
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    model = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    blocks = os.environ.get('BLOCKS', 'all')
    globals()['BLOCKS'] = BEHAV_BLOCKS if blocks == 'behav' else ALL_BLOCKS
    ckpt = f'results/fit_model{model}.json' if blocks == 'all' else f'results/fitb_model{model}.json'
    import calibrate as _c
    _c.BLOCKS = globals()['BLOCKS']
    E = build(); obj, tgt, scale = make_objective(E)
    P = json.load(open(ckpt)) if os.path.exists(ckpt) else dict(START)
    f = lambda v: obj(clip(unpack(v)), model=model)
    t0 = time.time()
    x0 = pack(P)
    step = float(os.environ.get('STEP', '0.12'))
    sim0 = np.vstack([x0] + [x0 + step * np.eye(len(x0))[i] for i in range(len(x0))])
    res = minimize(f, x0, method='Nelder-Mead',
                   options=dict(maxfev=steps, xatol=1e-4, fatol=1e-4, adaptive=True,
                                initial_simplex=sim0))
    P = clip(unpack(res.x))
    json.dump(P, open(ckpt, 'w'), indent=1)
    print(f"model {model}: objective {res.fun:.2f} after {res.nfev} evaluations "
          f"({(time.time()-t0)/60:.1f} min).  " + " ".join(f"{k}={P[k]:.4g}" for k in PARAM_NAMES))
