"""
targets.py -- the empirical quantities the models are screened and compared against.

Everything is read from the OSF archive.  Partial-activation counts are taken exactly
as published: myonset applies the 25 ms merge when it defines the activity windows, so
the event files already contain one event per merged window and must not be merged a
second time.

Burst durations are deliberately NOT used as targets: the integrated-profile step moves
onsets earlier and caps offsets at the threshold crossing, so a reported duration is the
excursion length plus an unknown onset displacement.
"""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from empirical import trials
from scipy import stats

CELLS = [(0.0, 0.3), (0.12, 0.3), (0.50, 0.3), (0.0, 2.0), (0.12, 2.0), (0.50, 2.0)]
LAB = ['0.3/c00', '0.3/c12', '0.3/c50', '2.0/c00', '2.0/c12', '2.0/c50']
BEH = os.path.join(os.environ.get('DELAY_DATA', '/home/claude/data'), 'data/behavioral')


def poisson_slope(y, x, iters=100):
    X = np.column_stack([np.ones_like(x), x])
    b = np.array([np.log(max(y.mean(), 1e-4)), 0.0])
    for _ in range(iters):
        mu = np.exp(np.clip(X @ b, -30, 30))
        z = X @ b + (y - mu) / np.maximum(mu, 1e-12)
        A = X.T @ (mu[:, None] * X) + 1e-9 * np.eye(2)
        bn = np.linalg.solve(A, X.T @ (mu * z))
        if np.max(np.abs(bn - b)) < 1e-11:
            return bn
        b = bn
    return b


def build(nsub=22):
    acc = np.zeros((nsub, 6)); pmet = np.zeros((nsub, 6)); met = np.zeros((nsub, 6))
    rt = np.zeros((nsub, 6)); rt_sd = np.zeros((nsub, 6))
    pmet_sd = np.zeros((nsub, 6)); met_sd = np.zeros((nsub, 6))
    dly = np.zeros((nsub, 6, 20)); stim = np.zeros((nsub, 6, 20)); ntr = np.zeros((nsub, 6))
    for s in range(nsub):
        A = {k: [] for k in range(6)}; P = {k: [] for k in range(6)}; M = {k: [] for k in range(6)}
        RT = {k: [] for k in range(6)}
        for tr in trials(s + 1):
            k = tr['dur'] * 3 + tr['coh']; ntr[s, k] += 1
            A[k].append(tr['acc']); RT[k].append(tr['rt'])
            P[k].append(tr['t_burst'] - tr['t_rs'])
            M[k].append(tr['t_resp'] - tr['t_burst'])
            v = tr['pma_t']
            w = v[(v >= tr['t_off'] + 1) & (v < tr['t_off'] + 3)] - tr['t_off'] - 1
            dly[s, k] += np.histogram(w, bins=20, range=(0, 2))[0]
            u = v[(v >= 0) & (v < 2.0) & (v < tr['t_burst'])]
            stim[s, k] += np.histogram(u, bins=20, range=(0, 2))[0]
        for k in range(6):
            acc[s, k] = np.mean(A[k]); pmet[s, k] = np.mean(P[k]); met[s, k] = np.mean(M[k])
            pmet_sd[s, k] = np.subtract(*np.percentile(P[k], [75, 25]))
            met_sd[s, k] = np.subtract(*np.percentile(M[k], [75, 25]))
            rt[s, k] = np.mean(RT[k])
            rt_sd[s, k] = np.subtract(*np.percentile(RT[k], [75, 25]))
    dly /= ntr[..., None]; stim /= ntr[..., None]

    # premature responses, from the behavioural files (response before the signal)
    prem = np.zeros((nsub, 6))
    for s in range(nsub):
        b = np.load(f'{BEH}/s{s+1}.npy', allow_pickle=True)
        coh = b[:, 5].astype(float); dur = b[:, 6].astype(float).astype(int)
        rm = b[:, 10].astype(float)
        for k, (c, D) in enumerate(CELLS):
            sel = (coh == c) & (dur == (18 if D == 0.3 else 120))
            prem[s, k] = (rm[sel] == 0).mean() if sel.any() else np.nan

    T = np.arange(20) * 0.1 + 0.05; T = T - T.mean()
    slope = np.full((nsub, 6), np.nan)
    for s in range(nsub):
        for k in range(6):
            if (dly[s, k] > 0).sum() >= 4:
                slope[s, k] = poisson_slope(dly[s, k], T)[1]

    lin, quad = delay_contrasts()

    return dict(acc=acc, rt=rt, rt_sd=rt_sd, pmet=pmet, met=met, pmet_sd=pmet_sd, met_sd=met_sd,
                dly=dly, stim=stim, prem=prem, slope=slope, trend=lin, quad=quad, ntr=ntr)


# --- the delay-period TREND, from the archive's own bin export -------------------------
# TREND (and QUAD) are read from analysis_PMA/data/PMA_rate_delay_JASP.csv rather than
# recomputed from `dly` above, for one reason: that export is the file the published
# analysis was run on, so contrasts taken from it reproduce the manuscript's F values
# exactly.  The event-file reconstruction in build() is a slightly different quantity --
# it gives a mean trend of 0.5551 against the export's 0.5604, a 1% difference coming from
# trial inclusion at the window edges -- and using it here would put the simulation's
# comparison target 1% away from the number the Results report.
JASP = os.path.join(os.environ.get('DELAY_DATA', '/home/claude/data'),
                    'analysis_PMA/data/PMA_rate_delay_JASP.csv')
BINW = 0.1                       # the export's bin width, seconds


def contrast_weights(x):
    """Weights of the linear and quadratic orthogonal-polynomial contrasts on the bin
    axis, normalised so that the contrast of a rate profile is a rate PER SECOND and per
    second squared -- i.e. the least-squares slope and curvature, not an arbitrary
    polynomial score.  Identical to fused_sweep.trend_contrasts, which is what the
    simulations report."""
    w1 = x / np.dot(x, x)
    q = x ** 2 - (x ** 2).mean()
    q = q - x * (np.dot(q, x) / np.dot(x, x))
    return w1, q / np.dot(q, q)


def delay_profiles(nsub=22):
    """(nsub, 6, 20) delay-period counts per trial per bin, read from the export.  Cell
    order is the LAB order: short duration at 0/12/50% coherence, then long."""
    import csv
    rows = list(csv.reader(open(JASP)))
    hdr = rows[0][1:]
    dat = np.array([[float(v) for v in r[1:]] for r in rows[1:]], float)
    cells = ['stim0_coh0', 'stim0_coh12', 'stim0_coh50',
             'stim1_coh0', 'stim1_coh12', 'stim1_coh50']
    X = np.stack([np.stack([dat[:, hdr.index(f'{c}_bin{b}')] for b in range(10, 30)], -1)
                  for c in cells], 1)
    return X[:nsub]


def delay_contrasts(nsub=22):
    """Per-participant, per-cell LINEAR (s^-2) and QUADRATIC (s^-3) contrasts of the
    delay-period rate.  Both are linear functionals of the rate, so -- unlike the
    log-rate slope -- the mean over cells and the contrast of the pooled profile agree."""
    try:
        X = delay_profiles(nsub)
    except (FileNotFoundError, OSError):
        return np.full((nsub, 6), np.nan), np.full((nsub, 6), np.nan)
    t = np.arange(20) * BINW + BINW / 2
    w1, w2 = contrast_weights(t - t.mean())
    r = X / BINW
    return r @ w1, r @ w2


def summary(E=None):
    E = E or build()
    return dict(
        acc=E['acc'].mean(0), rt=E['rt'].mean(0), rt_sd=E['rt_sd'].mean(0),
        pmet=E['pmet'].mean(0), met=E['met'].mean(0),
        met_sd=E['met_sd'].mean(0), pmet_sd=E['pmet_sd'].mean(0),
        dly=E['dly'].mean(0).mean(-1), stim=E['stim'].mean(0).mean(-1),
        prem=np.nanmean(E['prem'], 0), slope=np.nanmean(E['slope'], 0),
        dly_profile=E['dly'].mean(0))


if __name__ == '__main__':
    E = build(); S = summary(E)
    np.savez('targets.npz', **E)
    print(f"{'cell':>8} {'acc':>6} {'PMET':>6} {'sd':>6} {'MET':>6} {'sd':>6} "
          f"{'delay':>7} {'stim':>7} {'prem%':>6} {'slope':>6}")
    for k in range(6):
        print(f"{LAB[k]:>8} {S['acc'][k]:6.3f} {S['pmet'][k]:6.3f} {S['pmet_sd'][k]:6.3f} "
              f"{S['met'][k]:6.3f} {S['met_sd'][k]:6.3f} {S['dly'][k]:7.4f} "
              f"{S['stim'][k]:7.4f} {100*S['prem'][k]:6.2f} {S['slope'][k]:6.3f}")
