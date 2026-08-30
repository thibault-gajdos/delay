"""Behaviour-only simulation of the delayed-report model.

Identical dynamics to sim.simulate, but the gated-drive trace is not stored and the
window runs to the 4-s response deadline, so that the three response proportions the
calibration targets (pre-signal, faster than 0.2 s, no response by the deadline) are
all representable.  Returns per-trial motion direction, selected response and the
index at which the response criterion was first reached.
"""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim import DEFAULTS, DT

def simulate_behav(p, c, D, n, model=2, seed=0, deadline=4.25, gamma=0.0):
    q = dict(DEFAULTS); q.update(p)
    rng = np.random.default_rng(seed)
    lam, k, sx, sq = q['lam'], q['k'], q['sigx'], q['sigq']
    g, th, r = q['g'], q['theta'], q['r']
    kap, SDe = q['kappa'], q['SDe']
    se = SDe * np.sqrt(2 * kap)
    Lmot, LRS = q['Lmotion'], q['LRS']

    t_off = D + Lmot
    t_RS = D + q['delay'] + LRS
    Delta = t_RS - t_off
    beta = q['alpha'] * Delta if model == 2 else 0.0

    nt = int(round((t_RS + deadline) / DT)) + 1
    i_ev0, i_off, i_RS = (int(round(x / DT)) for x in (Lmot, t_off, t_RS))
    sqdt = np.sqrt(DT)

    s = rng.choice([-1.0, 1.0], size=n)
    x = np.zeros(n); y = np.zeros(n); eps = np.zeros(n); qq = np.zeros(n)
    m = np.zeros(n); sresp = np.zeros(n); maint = np.zeros(n, bool)
    yRS = np.zeros(n); resp_i = np.full(n, -1)

    for i in range(nt):
        t = i * DT
        ramp = 0.0 if i < i_off else (beta * (t - t_off) / Delta if i < i_RS else beta)
        yeff = y + sresp * ramp
        drive = np.where(sresp != 0, sresp * yeff - g, np.abs(yeff) - g)
        hit = (resp_i < 0) & (drive >= r)
        resp_i[hit] = i
        if i >= i_RS and np.all(resp_i >= 0):
            break
        if i < i_RS:
            ac = (~maint) & (i >= i_ev0) & (i < i_off)
            if ac.any():
                x[ac] += (s[ac] * k * c - gamma * x[ac]) * DT + sx * sqdt * rng.standard_normal(int(ac.sum()))
            inp = np.where(maint, m, x)
            if maint.any():
                eps[maint] += -kap * eps[maint] * DT + se * sqdt * rng.standard_normal(int(maint.sum()))
                inp = inp + np.where(maint, eps, 0.0)
            y = y + lam * (inp - y) * DT
            if i_ev0 <= i < i_off:
                cross = (~maint) & (np.abs(y) >= th)
                if cross.any():
                    maint[cross] = True; sresp[cross] = np.sign(y[cross])
                    m[cross] = sresp[cross] * th; y[cross] = sresp[cross] * th
            if i == i_off - 1:
                left = ~maint
                if left.any():
                    sg = np.sign(y[left]); sg[sg == 0] = 1.0
                    maint[left] = True; sresp[left] = sg; m[left] = y[left]
            if i == i_RS - 1:
                yRS = y.copy()
        else:
            qq += (q['muRS'] - gamma * qq) * DT + sq * sqdt * rng.standard_normal(n)
            y = y + lam * (yRS + sresp * qq - y) * DT
    return dict(s=s, sresp=sresp, m=m, resp_i=resp_i,
                i_phys=int(round((D + q['delay']) / DT)), Lmotor=q['Lmotor'])


def behav_rt_split(out):
    """Correct- and incorrect-response RT arrays on valid trials, plus the premature
    count, for the quantile+premature-rate calibration (2026-08-30: premature responses
    are scored again, as a proportion alongside correct/incorrect -- see
    calib_behav_q.py's docstring). <0.2 s and no-response-by-deadline trials are still
    dropped entirely, not scored, not counted in any denominator: the historical finding
    that no-response-by-deadline is not reproducible at any behaviourally-plausible
    setting (README, "Removed 2026-08-30") is untouched by this change."""
    ri, i_phys, Lm = out['resp_i'], out['i_phys'], out['Lmotor']
    responded = ri >= 0
    prem = responded & (ri < i_phys)
    rt = np.where(responded & ~prem, (ri - i_phys) * DT + Lm, np.nan)
    valid = responded & ~prem & (rt >= 0.2) & (rt <= 4.0)
    correct = valid & (out['s'] == out['sresp'])
    incorrect = valid & (out['s'] != out['sresp'])
    return rt[correct], rt[incorrect], int(prem.sum())

