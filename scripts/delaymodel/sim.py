"""
sim.py -- the extended evidence-accumulation model for delayed reports, and the rule that
turns the simulated drive into partial activations.

No EMG post-processing is applied to the simulation.  An activation begins at every
upcrossing of the gate by the drive of the selected channel: no merging of nearby
excursions, no minimum duration, no refinement of onsets.  The excursion that reaches the
response criterion r is the response burst; its onset is where premotor time is measured
from, and it is not counted as a partial activation.

Scale convention: sigma_x fixes the amplitude unit; lambda fixes the time unit of the
filter.  All other quantities are free.
"""
import numpy as np

DT = 1.0 / 1024.0    # the EMG/EEG acquisition rate

# Maintenance-noise initialisation.  WARM_START=True draws eps from its stationary law
# N(0, SDe^2) at the moment a trial enters maintenance; False leaves eps at 0.  H(tau)=0
# either way -- the filter output is continuous and m is defined as its value at tau --
# so only eps is at issue.  Starting eps at zero is the sole source of a slow (1/kappa)
# transient in the model; started at stationarity the moments of H relax at rate at
# least lambda, uniformly in kappa, so the delay-period analysis window lies entirely
# outside the transient (Appendix, Proposition A.1).  Kept switchable for that
# comparison; the manuscript uses WARM_START=True throughout.
WARM_START = True

PARAM_NAMES = ['k', 'lam', 'theta', 'g', 'kappa', 'SDe', 'alpha', 'muRS',
               'sigq', 'r', 'Lmotor', 'Lmotion', 'LRS']

DEFAULTS = dict(sigx=0.100, delay=3.0)


def excursions(active):
    """(start, stop) indices of the maximal runs of True in a boolean array: the
    excursions of the drive above the gate, each one activation."""
    d = np.diff(np.concatenate([[0], active.astype(np.int8), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def simulate(p, c, D, n, model=2, seed=0, post=1.8, gamma_x=0.0, gamma_q=0.0):
    """One condition.  `p` is a dict of the free parameters (PARAM_NAMES).

    Returns a dict with, per trial: the selected side, the maintained level, the
    response-trigger index, and the gated drive of the selected channel.
    `gamma_x` is the leak of the motion-evidence accumulator x, `gamma_q` the leak of
    the response-signal accumulator q; both are zero in the model as specified."""
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

    nt = int(round((t_RS + post) / DT)) + 1
    i_ev0, i_off, i_RS = (int(round(x / DT)) for x in (Lmot, t_off, t_RS))
    sqdt = np.sqrt(DT)

    s = rng.choice([-1.0, 1.0], size=n)
    x = np.zeros(n); y = np.zeros(n); eps = np.zeros(n); qq = np.zeros(n)
    m = np.zeros(n); sresp = np.zeros(n); maint = np.zeros(n, bool)
    yRS = np.zeros(n); resp_i = np.full(n, -1)
    Z = np.empty((n, nt), np.float32)          # gated drive of the selected channel

    for i in range(nt):
        t = i * DT
        if i < i_off:
            ramp = 0.0
        elif i < i_RS:
            ramp = beta * (t - t_off) / Delta
        else:
            ramp = beta
        yeff = y + sresp * ramp
        drive = np.where(sresp != 0, sresp * yeff - g, np.abs(yeff) - g)
        Z[:, i] = drive
        hit = (resp_i < 0) & (drive >= r)
        resp_i[hit] = i
        if i >= i_RS and np.all(resp_i >= 0):
            Z = Z[:, :i + 1]; nt = i + 1
            break

        if i < i_RS:
            acc = (~maint) & (i >= i_ev0) & (i < i_off)
            if acc.any():
                na = int(acc.sum())
                x[acc] += (s[acc] * k * c - gamma_x * x[acc]) * DT + sx * sqdt * rng.standard_normal(na)
            inp = np.where(maint, m, x)
            if maint.any():
                nm = int(maint.sum())
                eps[maint] += -kap * eps[maint] * DT + se * sqdt * rng.standard_normal(nm)
                inp = inp + np.where(maint, eps, 0.0)
            y = y + lam * (inp - y) * DT
            if i_ev0 <= i < i_off:
                cross = (~maint) & (np.abs(y) >= th)
                if cross.any():
                    maint[cross] = True
                    sresp[cross] = np.sign(y[cross])
                    m[cross] = sresp[cross] * th
                    y[cross] = sresp[cross] * th
                    if WARM_START:
                        eps[cross] = SDe * rng.standard_normal(int(cross.sum()))
            if i == i_off - 1:
                left = ~maint
                if left.any():
                    sg = np.sign(y[left]); sg[sg == 0] = 1.0
                    maint[left] = True; sresp[left] = sg; m[left] = y[left]
                    if WARM_START:
                        eps[left] = SDe * rng.standard_normal(int(left.sum()))
            if i == i_RS - 1:
                yRS = y.copy()
        else:
            qq += (q['muRS'] - gamma_q * qq) * DT + sq * sqdt * rng.standard_normal(n)
            y = y + lam * (yRS + sresp * qq - y) * DT

    return dict(Z=Z, s=s, sresp=sresp, m=m, resp_i=resp_i, i_RS=i_RS, i_off=i_off,
                D=D, t_RS=t_RS, nt=nt, p=q)


def measures(sim, windows=True):
    """Behaviour, latencies and partial-activation onset times, relative to stimulus
    onset.  An activation begins at each upcrossing of the gate (drive > 0); the last one
    at or before the response index is the onset of the response burst."""
    Z, ri, sr, s = sim['Z'], sim['resp_i'], sim['sresp'], sim['s']
    q = sim['p']; D, i_RS = sim['D'], sim['i_RS']
    i_phys = int(round((D + q['delay']) / DT))          # physical RS onset
    n, nt = Z.shape
    rt = np.where(ri >= 0, (ri - i_phys) * DT + q['Lmotor'], np.nan)
    premature = (ri >= 0) & (ri < i_phys)
    noresp = ri < 0
    valid = (~premature) & (~noresp) & (rt >= 0.2) & (rt <= 4.0)

    act = Z > 0
    up = np.empty_like(act)
    up[:, 0] = act[:, 0]
    up[:, 1:] = act[:, 1:] & ~act[:, :-1]
    stop = np.where(ri >= 0, ri, nt - 1)
    up &= np.arange(nt)[None, :] <= stop[:, None]
    rows, cols = np.nonzero(up)                        # row-major, so sorted by trial
    bounds = np.searchsorted(rows, np.arange(n + 1))

    pmet = np.full(n, np.nan); met = np.full(n, np.nan)
    events = []
    for j in range(n):
        c = cols[bounds[j]:bounds[j + 1]]
        if ri[j] >= 0 and c.size:
            onset = c[-1]
            pmet[j] = (onset - i_phys) * DT
            met[j] = (ri[j] - onset) * DT + q['Lmotor']
            c = c[:-1]
        events.append(c * DT if valid[j] else np.empty(0))
    return dict(acc=(s == sr), rt=rt, pmet=pmet, met=met, valid=valid,
                premature=premature, noresp=noresp, events=events, D=D)


def rate_profile(mm, t0, t1, nbin):
    """Mean number of events per bin per valid trial, window [t0,t1) from onset."""
    v = np.flatnonzero(mm['valid'])
    if v.size == 0:
        return np.zeros(nbin)
    h = np.zeros(nbin)
    for j in v:
        e = mm['events'][j]
        h += np.histogram(e[(e >= t0) & (e < t1)], bins=nbin, range=(t0, t1))[0]
    return h / v.size
