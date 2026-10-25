"""
targets.py -- the empirical quantities the models are compared against, per participant
and per coherence x duration condition, read from the OSF archive (DELAY_DATA).

Accuracy, response time, premotor time (PMT) and motor time (MT), and the PMA counts of the
stimulus-onset and delay windows, are taken from the trial-level event files
(empirical.trials).  The slope and curvature of the PMA rate in the delay window, and its
slope in the stimulus-onset window, are taken from the archive's own bin exports
(analysis_PMA/data/PMA_rate_*_JASP.csv), the files the statistical analysis runs on, so that
the observed values the models are compared against are those of the Results.

PMA counts are used as exported: the detection procedure merges intervals separated by
less than 25 ms when it defines the activity windows, so the event files contain one event
per merged window.
"""
import numpy as np, os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from empirical import trials

CELLS = [(0.0, 0.3), (0.12, 0.3), (0.50, 0.3), (0.0, 2.0), (0.12, 2.0), (0.50, 2.0)]
JASPDIR = os.path.join(os.environ.get('DELAY_DATA', '/home/claude/data'), 'analysis_PMA/data')
BINW = 0.1                       # the exports' bin width, seconds
# the exported windows: file, and the bin indices it uses
WINDOWS = dict(stim=('PMA_rate_stim_JASP.csv', range(0, 20)),
               delay=('PMA_rate_delay_JASP.csv', range(10, 30)))


def build(nsub=22):
    """Per participant (rows) and condition (columns, order of CELLS): accuracy, mean RT,
    mean PMT and MT in s; PMA counts per trial in the 20 bins of the stimulus-onset and delay
    windows; and the slopes (s^-2) and curvature (s^-3) of the exported PMA rates."""
    acc = np.zeros((nsub, 6)); rt = np.zeros((nsub, 6))
    pmet = np.zeros((nsub, 6)); met = np.zeros((nsub, 6))
    dly = np.zeros((nsub, 6, 20)); stim = np.zeros((nsub, 6, 20)); ntr = np.zeros((nsub, 6))
    for s in range(nsub):
        A = {k: [] for k in range(6)}; RT = {k: [] for k in range(6)}
        P = {k: [] for k in range(6)}; M = {k: [] for k in range(6)}
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
            acc[s, k] = np.mean(A[k]); rt[s, k] = np.mean(RT[k])
            pmet[s, k] = np.mean(P[k]); met[s, k] = np.mean(M[k])
    dly /= ntr[..., None]; stim /= ntr[..., None]
    trend, quad = window_contrasts('delay', nsub)
    stim_trend, _ = window_contrasts('stim', nsub)
    return dict(acc=acc, rt=rt, pmet=pmet, met=met, dly=dly, stim=stim, ntr=ntr,
                trend=trend, quad=quad, stim_trend=stim_trend)


def contrast_weights(x):
    """Weights giving the linear and quadratic coefficients of the least-squares fit
    y = b0 + b1 x + b2 x^2 on centred times x (seconds); identical to
    fused_sweep.trend_contrasts, which the simulations use."""
    w1 = x / np.dot(x, x)
    q = x ** 2 - (x ** 2).mean()
    q = q - x * (np.dot(q, x) / np.dot(x, x))
    return w1, q / np.dot(q, q)


def window_profiles(window, nsub=22):
    """(nsub, 6, nbin) PMA counts per trial per bin of one window, from the export."""
    fname, bins = WINDOWS[window]
    rows = list(csv.reader(open(os.path.join(JASPDIR, fname))))
    hdr = rows[0][1:]
    dat = np.array([[float(v) for v in r[1:]] for r in rows[1:]], float)
    cells = ['stim0_coh0', 'stim0_coh12', 'stim0_coh50',
             'stim1_coh0', 'stim1_coh12', 'stim1_coh50']
    X = np.stack([np.stack([dat[:, hdr.index(f'{c}_bin{b}')] for b in bins], -1)
                  for c in cells], 1)
    return X[:nsub]


def window_contrasts(window, nsub=22):
    """Per-participant, per-condition slope (s^-2) and curvature (s^-3) of the PMA rate
    (counts per trial divided by the bin width) of one exported window."""
    X = window_profiles(window, nsub)
    n = X.shape[-1]
    t = np.arange(n) * BINW + BINW / 2
    w1, w2 = contrast_weights(t - t.mean())
    r = X / BINW
    return r @ w1, r @ w2
