"""
figures.py -- the figures that show simulated dynamics.

Nothing is fitted, so these are drawn at the REPRESENTATIVE ADMISSIBLE DRAW selected by
representative.py (see that file for the rule and for why the median of the ranges is not
usable). The same parameter vector is used for both models, so the two trajectory figures
differ only in beta -- which is the only term in which the architectures differ.

Run representative.py first.  Usage: python3 figures.py
"""
import numpy as np, json, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, rate_profile, DT, _bursts
from delaymodel.targets import CELLS, LAB, build, poisson_slope

COL = {0: '#9ecae1', 12: '#3182bd', 50: '#08306b'}
CC = [0, 12, 50]
STY = {0.3: '-', 2.0: '--'}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'lines.linewidth': 1.1})


def cell_col(k):
    return COL[CC[k % 3]]


def cell_sty(k):
    return '-' if k < 3 else '--'


# ---------------------------------------------------------------- trajectories
def fig_trajectories(P, model=2, n=800, ntraj=18, fname='plots/fig_traj.pdf'):
    """Simulated counterpart of the manuscript's Figure 1: the gated motor-preparatory
    signal, time-locked to the four trial events.

    Individual trajectories rather than a standard-deviation band. A band is a cross-section
    -- it says how far apart trials are at each instant, which is not what one wants to know
    about a process whose whole content is what a trial DOES over time. Twenty single trials
    drawn in grey behind the six condition means show the same variability in the form that
    carries the mechanism: the excursions towards the gate, their duration, and the fact
    that a trial which crosses does so as an event rather than as a displacement of a mean.
    """
    fig, ax = plt.subplots(1, 4, figsize=(9.2, 2.7),
                           gridspec_kw=dict(width_ratios=[2.2, 3, 0.9, 0.9], wspace=0.12))
    # --- individual trials, pooled over conditions, drawn first and behind everything
    st = simulate(P, 0.12, 2.0, ntraj * 3, model, seed=17)
    ms = measures(st)
    Yt = st['Z'] + P['g']; rit = st['resp_i']
    tt = np.arange(Yt.shape[1]) * DT
    picked = [j for j in np.flatnonzero(ms['valid'])][:ntraj]
    for j in picked:
        y = np.where(np.arange(Yt.shape[1]) <= rit[j], Yt[j], np.nan)
        ax[0].plot(tt[tt <= 2.2], y[tt <= 2.2], color='0.55', lw=.35, alpha=.55, zorder=0)
        t2 = tt - 2.0; w2 = (t2 >= -0.1) & (t2 <= 3 + P['LRS'])
        ax[1].plot(t2[w2], y[w2], color='0.55', lw=.35, alpha=.55, zorder=0)
        t3 = tt - 5.0; w3 = (t3 >= 0) & (t3 <= 0.45)
        ax[2].plot(t3[w3], y[w3], color='0.55', lw=.35, alpha=.55, zorder=0)
        L = int(0.25 / DT)
        if rit[j] >= L:
            ax[3].plot(np.arange(-L, 0) * DT, Yt[j, rit[j] - L:rit[j]],
                       color='0.55', lw=.35, alpha=.55, zorder=0)
    # --- condition means
    for k, (c, D) in enumerate(CELLS):
        s_ = simulate(P, c, D, n, model, seed=41 + k)
        m = measures(s_)
        Y = s_['Z'] + P['g']
        v = m['valid']; ri = s_['resp_i']
        col, sty = cell_col(k), cell_sty(k)
        t = np.arange(Y.shape[1]) * DT
        M = np.array([np.where(np.arange(Y.shape[1]) <= (ri[j] if ri[j] >= 0 else 10**9),
                               Y[j], np.nan) for j in np.flatnonzero(v)])
        mu = np.nanmean(M, 0)
        w = t <= 2.2
        ax[0].plot(t[w], mu[w], sty, color=col, lw=1.4, zorder=3)
        t2 = t - D
        w2 = (t2 >= -0.1) & (t2 <= 3 + P['LRS'])
        ax[1].plot(t2[w2], mu[w2], sty, color=col, lw=1.4, zorder=3)
        t3 = t - (D + 3)
        w3 = (t3 >= 0) & (t3 <= 0.45)
        ax[2].plot(t3[w3], mu[w3], sty, color=col, lw=1.4, zorder=3)
        L = int(0.25 / DT)
        R = np.full((v.sum(), L), np.nan)
        for i, j in enumerate(np.flatnonzero(v)):
            if ri[j] >= L:
                R[i] = Y[j, ri[j] - L:ri[j]]
        ax[3].plot(np.arange(-L, 0) * DT, np.nanmean(R, 0), sty, color=col, lw=1.4, zorder=3)
    for a, ttl in zip(ax, ['from stimulus onset', 'from stimulus offset',
                           'from response signal', 'from response']):
        a.axhline(P['g'], color='0.3', ls=':', lw=.8)
        a.axhline(P['theta'], color='0.6', ls=':', lw=.8)
        a.set_title(ttl, fontsize=8)
        a.set_ylim(-0.01, P['g'] + P['r'] * 1.6)
    ax[0].set_ylabel('gated motor-preparatory signal')
    ax[0].text(0.02, P['g'] + 0.004, '$g$', fontsize=8, color='0.3')
    ax[0].text(0.02, P['theta'] + 0.004, r'$\theta$', fontsize=8, color='0.6')
    for a in ax[1:]:
        a.set_yticklabels([])
    h = [plt.Line2D([], [], color=COL[c], ls='-') for c in CC] + \
        [plt.Line2D([], [], color='0.4', ls=s) for s in ['-', '--']]
    ax[1].legend(h, ['0%', '12%', '50%', '0.3 s', '2 s'], ncol=5, fontsize=6.5,
                 frameon=False, loc='upper left')
    ax[0].legend([plt.Line2D([], [], color='0.55', lw=.5)], ['single trials'],
                 fontsize=6.5, frameon=False, loc='upper left')
    ax[0].set_xlabel('from onset (s)'); ax[1].set_xlabel('from offset (s)')
    ax[2].set_xlabel('(s)'); ax[3].set_xlabel('(s)')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- single trial
def fig_single(P, n=400, fname='plots/fig_single.pdf'):
    s = simulate(P, 0.12, 2.0, n, 2, seed=7)
    m = measures(s)
    cand = [j for j in np.flatnonzero(m['valid'])
            if len(m['events'][j]) >= 1 and any(e > s['D'] + 1 for e in m['events'][j])]
    j = cand[0]
    Y = s['Z'][j] + P['g']; t = np.arange(len(Y)) * DT
    ri = s['resp_i'][j]
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 3.2), sharex=True,
                           gridspec_kw=dict(height_ratios=[3, 1], hspace=.12))
    ax[0].plot(t[:ri + 1], Y[:ri + 1], color='#08306b', lw=.7)
    ax[0].axhline(P['g'], color='0.3', ls=':', lw=.8)
    ax[0].axhline(P['g'] + P['r'], color='#c1272d', ls=':', lw=.8)
    ax[0].text(0.05, P['g'] + 0.002, '$g$', fontsize=8, color='0.3')
    ax[0].text(0.05, P['g'] + P['r'] + 0.002, '$g+r$', fontsize=8, color='#c1272d')
    for x, lab in [(s['D'], 'offset'), (s['D'] + 3, 'signal')]:
        for a in ax:
            a.axvline(x, color='0.7', lw=.7)
        ax[0].text(x + .02, ax[0].get_ylim()[1] * .95, lab, fontsize=6.5, color='0.5')
    act = Y[:ri + 1] > P['g']
    for a, b in _bursts(act):
        col = '#c1272d' if b >= ri - 1 else '#e6892f'
        ax[1].axvspan(a * DT, b * DT, color=col, alpha=.85)
    ax[1].set_yticks([]); ax[1].set_ylabel('detected', fontsize=7)
    ax[1].set_xlabel('time from stimulus onset (s)')
    ax[0].set_ylabel('gated signal')
    # panel letter: this figure is panel B of the manuscript's Figure 2, whose panel A is a
    # hand-drawn schematic; the letter used to be added by hand after export
    fig.text(0.005, 0.94, 'B', fontsize=20, weight='bold', va='top')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- rates
def fig_rates(P, E, n=900, fname='plots/fig_rates.pdf'):
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 2.7), sharey=True,
                           gridspec_kw=dict(width_ratios=[1, 1], wspace=.08))
    for k, (c, D) in enumerate(CELLS):
        s = simulate(P, c, D, n, 2, seed=131 + k); m = measures(s)
        md = rate_profile(m, D + 1, D + 3, 20)
        ms = rate_profile(m, 0.0, 2.0, 20)
        tb = np.arange(20) * 0.1 + 0.05
        col, sty = cell_col(k), cell_sty(k)
        ax[0].plot(tb, E['stim'].mean(0)[k], sty, color=col, alpha=.35, lw=.8)
        ax[0].plot(tb, ms, sty, color=col, lw=1.4)
        ax[1].plot(tb + 1, E['dly'].mean(0)[k], sty, color=col, alpha=.35, lw=.8)
        ax[1].plot(tb + 1, md, sty, color=col, lw=1.4)
    ax[0].set_title('stimulus period', fontsize=8)
    ax[1].set_title('delay period', fontsize=8)
    ax[0].set_ylabel('partial activations per 0.1 s bin')
    ax[0].set_xlabel('time from stimulus onset (s)')
    ax[1].set_xlabel('time from stimulus offset (s)')
    h = [plt.Line2D([], [], color='0.4', alpha=.35, lw=.8), plt.Line2D([], [], color='0.4', lw=1.4)]
    ax[1].legend(h, ['data', 'model'], fontsize=6.5, frameon=False)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- shape tests
def fig_shape(E, fname='plots/fig_shape.pdf'):
    tb = np.arange(20) * 0.1 + 0.05; tc = tb - tb.mean()
    lr = np.log(E['dly'].mean(0) .mean(-1) / 0.1)
    sl = np.nanmean(E['slope'], 0)
    fig, ax = plt.subplots(1, 3, figsize=(9.0, 2.5), gridspec_kw=dict(wspace=.30))
    for k in range(6):
        y = E['dly'].mean(0)[k]
        ax[0].plot(tb + 1, np.log(y / 0.1), cell_sty(k), color=cell_col(k))
    ax[0].set_xlabel('time from stimulus offset (s)'); ax[0].set_ylabel(r'$\ln\nu$')
    ax[0].set_title('delay log-rate profiles', fontsize=8)
    for k in range(6):
        ax[1].plot(lr[k], sl[k], 'o', color=cell_col(k),
                   mfc=cell_col(k) if k >= 3 else 'w', ms=5)
    b = np.polyfit(lr, sl, 1)
    xs = np.linspace(lr.min(), lr.max(), 10)
    ax[1].plot(xs, np.polyval(b, xs), '-', color='0.5', lw=.8)
    ax[1].set_xlabel(r'$\ln\nu$'); ax[1].set_ylabel('log-rate slope (s$^{-1}$)')
    ax[1].set_title(f'P2: slope falls as rate rises', fontsize=8)
    for k in range(6):
        ax[2].plot(sl[k] ** 2, lr[k], 'o', color=cell_col(k),
                   mfc=cell_col(k) if k >= 3 else 'w', ms=5)
    b2 = np.polyfit(sl ** 2, lr, 1)
    xs = np.linspace((sl ** 2).min(), (sl ** 2).max(), 10)
    ax[2].plot(xs, np.polyval(b2, xs), '-', color='0.5', lw=.8)
    ax[2].set_xlabel('(log-rate slope)$^2$'); ax[2].set_ylabel(r'$\ln\nu$')
    ax[2].set_title(r'curvature-free identity, slope $=-A/2b^2$', fontsize=8)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# --------------------------------------------------- out-of-sample delay slope
# Cell labels are duration/coherence ("0.3/c12" = 0.3 s duration, 12% coherence);
# the out-of-sample figures group by duration and spell coherence out instead, since
# that compact form reads poorly on an axis.
DUR_LAB = {0.3: '0.3 s', 2.0: '2.0 s'}
COH_LAB = {0.0: '0%', 0.12: '12%', 0.50: '50%'}
# ---------------------------------------------------------------- Model 1
if __name__ == '__main__':
    R = json.load(open('results/representative.json'))
    P1 = dict(R['params'], alpha=0.0)
    P2 = dict(R['params'], alpha=R['alpha'])
    E = build()
    fig_trajectories(P2, 2, fname='plots/fig_traj_m2.pdf')
    fig_trajectories(P1, 1, fname='plots/fig_traj_m1.pdf')
    fig_single(P2)
    fig_rates(P2, E)
    fig_shape(E)
    print(f"figures written at the representative draw {R['draw']} "
          f"(beta = {R['beta']:.4f}, the Model 2 admissible median)")
