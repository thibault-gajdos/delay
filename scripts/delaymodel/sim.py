"""
sim.py -- the extended evidence-accumulation model for delayed reports, exactly as
specified in Achard et al. (R1 v3), together with the detection rules that turn the
simulated gated drive into countable partial activations.

Detection follows myonset (Spieser), detector_var, as used for the data:
    th = 7 SD above baseline          -> here: the gated drive exceeds zero
    time_limit = 0.025 s              -> intervals separated by <= 25 ms are merged
    min_samples = 3, varying_min = 0  -> merged intervals with < 3 supra-threshold
                                         samples are discarded
The integrated-profile refinement is not simulated: it moves the boundaries of each
interval but does not change their number, and only the number is counted.

Scale convention: sigma_x fixes the amplitude unit; lambda fixes the time unit of the
filter.  All other quantities are free.
"""
import numpy as np

DT = 1.0 / 1024.0    # native EMG/EEG acquisition rate (was 1 ms; see myonset below)
MERGE_S = round(0.025 / DT)   # myonset time_limit = 0.025 s, in samples at DT
MIN_SAMP = 3          # myonset min_samples, already a native-rate sample count

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


def _bursts(active, merge=MERGE_S, min_samp=MIN_SAMP, varying=True):
    """Indices (start, stop) of merged supra-threshold intervals of a boolean array.

    Merges gaps <= `merge`, then drops intervals with fewer than the minimum number of
    supra-threshold samples counted inside the merged interval.  With `varying` the
    minimum is adaptive, as in myonset with varying_minRaw = 1: it is
    `min_samp + f`, where f is the number of intervals per second of epoch containing
    at most `2*min_samp` supra-threshold samples."""
    if not active.any():
        return []
    d = np.diff(active.astype(np.int8))
    st = list(np.flatnonzero(d == 1) + 1)
    en = list(np.flatnonzero(d == -1) + 1)
    if active[0]:
        st = [0] + st
    if active[-1]:
        en = en + [len(active)]
    out = [[st[0], en[0]]]
    for a, b in zip(st[1:], en[1:]):
        if a - out[-1][1] <= merge:
            out[-1][1] = b
        else:
            out.append([a, b])
    cnt = np.array([active[a:b].sum() for a, b in out])
    if varying:
        f = (cnt <= 2 * min_samp).sum() / (len(active) * DT)
        thr = min_samp + f
    else:
        thr = min_samp
    return [(a, b) for (a, b), c in zip(out, cnt) if c >= thr]


def simulate(p, c, D, n, model=2, seed=0, post=1.8, gamma=0.0):
    """One condition.  `p` is a dict of the free parameters (PARAM_NAMES).

    Returns a dict with, per trial: the selected side, the maintained level, the
    response-trigger index, and the gated drive of the selected channel.
    `gamma` adds a common leak to both accumulators (reviewer R2.3)."""
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
                x[acc] += (s[acc] * k * c - gamma * x[acc]) * DT + sx * sqdt * rng.standard_normal(na)
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
            qq += (q['muRS'] - gamma * qq) * DT + sq * sqdt * rng.standard_normal(n)
            y = y + lam * (yRS + sresp * qq - y) * DT

    return dict(Z=Z, s=s, sresp=sresp, m=m, resp_i=resp_i, i_RS=i_RS, i_off=i_off,
                D=D, t_RS=t_RS, nt=nt, p=q)


def measures(sim, windows=True):
    """Behaviour, latencies and partial-activation event times, with the detection
    rules applied.  Times are returned relative to stimulus onset."""
    Z, ri, sr, s = sim['Z'], sim['resp_i'], sim['sresp'], sim['s']
    q = sim['p']; D, i_RS, nt = sim['D'], sim['i_RS'], sim['nt']
    i_phys = int(round((D + q['delay']) / DT))          # physical RS onset
    n = Z.shape[0]
    rt = np.where(ri >= 0, (ri - i_phys) * DT + q['Lmotor'], np.nan)
    premature = (ri >= 0) & (ri < i_phys)
    noresp = ri < 0
    valid = (~premature) & (~noresp) & (rt >= 0.2) & (rt <= 4.0)

    pmet = np.full(n, np.nan); met = np.full(n, np.nan)
    events = []
    for j in range(n):
        stop = ri[j] + 1 if ri[j] >= 0 else nt
        act = Z[j, :stop] > 0
        bs = _bursts(act)
        if ri[j] >= 0 and bs:
            onset = bs[-1][0]
            pmet[j] = (onset - i_phys) * DT
            met[j] = (ri[j] - onset) * DT + q['Lmotor']
            bs = bs[:-1]
        if valid[j]:
            events.append(np.array([b[0] for b in bs]) * DT)
        else:
            events.append(np.empty(0))
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
