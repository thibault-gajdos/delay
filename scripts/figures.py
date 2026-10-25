"""
figures.py -- the figures that show simulated dynamics: manuscript Figure 1 (the mean
signal reaching the gate under each model, over thirty parameter vectors admissible in both
models; ../manuscript/latex/figures/fig01.png and fig02.png) and Figure 2B (one simulated
trial at the vector chosen by representative.py; fig04.png).

Run representative.py first.  Usage: python3 figures.py
"""
import numpy as np, json, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT, excursions
from delaymodel.targets import CELLS

COL = {0: '#9ecae1', 12: '#3182bd', 50: '#08306b'}
C_M1, C_M2 = '#c1272d', '#3182bd'
# common time grids for Figure 1, so every vector's mean can be interpolated onto the same axis
TG0 = np.linspace(0.0, 2.2, 300)
TG1 = np.linspace(-0.1, 3.2, 400)
TG2 = np.linspace(0.0, 0.45, 120)
CC = [0, 12, 50]
STY = {0.3: '-', 2.0: '--'}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'lines.linewidth': 1.1})


def cell_col(k):
    return COL[CC[k % 3]]


def cell_sty(k):
    return '-' if k < 3 else '--'


# ---------------------------------------------------------------- trajectories
def fig_trajectories(model=2, ndraw=30, n=300, fname='plots/fig_traj.pdf', png=None,
                     dump='results/sweep_screened.jsonl'):
    """Manuscript Figure 1: the signal reaching the gate, time-locked to the four trial
    events. Each thin line is the mean over trials of one admissible parameter vector; the
    thick line is the mean over vectors. The vertical axis is in units of each vector's own
    distance from the maintenance target to the gate (0: maintenance target, 1: gate g).
    The vectors are admissible in both models and are the same in both panels, so the two
    panels differ only in beta.
    """
    recs = [json.loads(l) for l in open(dump)]
    both = [r for r in recs if r['model1']['admissible'] and r['model2']['admissible']]
    step = max(1, len(both) // ndraw)
    sel = both[::step][:ndraw]
    fig, ax = plt.subplots(1, 4, figsize=(9.2, 2.7),
                           gridspec_kw=dict(width_ratios=[2.2, 3, 0.9, 0.9], wspace=0.12))
    L = int(0.25 / DT)
    acc = {0: [], 1: [], 2: [], 3: []}
    for q, rec in enumerate(sel):
        P = dict(rec['params'], alpha=(rec['alpha'] if model == 2 else 0.0))
        curves = {0: [], 1: [], 2: [], 3: []}
        for k, (c, D) in enumerate(CELLS):
            s_ = simulate(P, c, D, n, model, seed=8000 + 97 * q + k)
            m = measures(s_)
            v = np.flatnonzero(m['valid'])
            if v.size == 0:
                continue
            Y = s_['Z'] + P['g']; ri = s_['resp_i']
            mbar = np.nanmean(np.abs(s_['m'][v]))
            scale = max(P['g'] - mbar, 1e-6)
            M = np.array([np.where(np.arange(Y.shape[1]) <= ri[j], Y[j], np.nan) for j in v])
            mu = (np.nanmean(M, 0) - mbar) / scale
            t = np.arange(Y.shape[1]) * DT
            curves[0].append(np.interp(TG0, t, mu, left=np.nan, right=np.nan))
            curves[1].append(np.interp(TG1, t - D, mu, left=np.nan, right=np.nan))
            curves[2].append(np.interp(TG2, t - (D + 3), mu, left=np.nan, right=np.nan))
            R = np.full((v.size, L), np.nan)
            for a_, j in enumerate(v):
                if ri[j] >= L:
                    R[a_] = (Y[j, ri[j] - L:ri[j]] - mbar) / scale
            curves[3].append(np.nanmean(R, 0))
        col = C_M1 if model == 1 else C_M2
        for w, g in enumerate([TG0, TG1, TG2, np.arange(-L, 0) * DT]):
            if not curves[w]:
                continue
            y = np.nanmean(np.array(curves[w]), 0)
            acc[w].append(y)
            ax[w].plot(g, y, color=col, lw=.5, alpha=.45, zorder=1)
    for w, g in enumerate([TG0, TG1, TG2, np.arange(-L, 0) * DT]):
        if acc[w]:
            ax[w].plot(g, np.nanmean(np.array(acc[w]), 0),
                       color=(C_M1 if model == 1 else C_M2), lw=2.0, zorder=3)
    for a, ttl in zip(ax, ['from stimulus onset', 'from stimulus offset',
                           'from response signal', 'from response']):
        a.axhline(1.0, color='0.3', ls=':', lw=.8)
        a.axhline(0.0, color='0.6', ls=':', lw=.8)
        a.set_title(ttl, fontsize=8)
        a.set_ylim(-0.6, 1.9)
    ax[0].set_ylabel('signal, in units of the distance\nfrom the maintenance target to the gate')
    ax[0].text(0.03, 1.04, '$g$', fontsize=8, color='0.3')
    ax[0].text(0.03, 0.04, 'maintenance target', fontsize=7, color='0.55')
    for a in ax[1:]:
        a.set_yticklabels([])
    ax[1].legend([plt.Line2D([], [], color=(C_M1 if model == 1 else C_M2), lw=.5),
                  plt.Line2D([], [], color=(C_M1 if model == 1 else C_M2), lw=2.0)],
                 [f'one admissible parameter vector ({len(sel)} shown)', 'mean over vectors'],
                 fontsize=6.5, frameon=False, loc='upper left')
    ax[0].set_xlabel('from onset (s)'); ax[1].set_xlabel('from offset (s)')
    ax[2].set_xlabel('(s)'); ax[3].set_xlabel('(s)')
    fig.savefig(fname, bbox_inches='tight')
    if png:
        fig.savefig(png, dpi=200, bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- single trial
def fig_single(P, n=400, fname='plots/fig_single.pdf', png=None):
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
    for a, b in excursions(act):
        col = '#c1272d' if b >= ri - 1 else '#e6892f'
        ax[1].axvspan(a * DT, b * DT, color=col, alpha=.85)
    ax[1].set_yticks([]); ax[1].set_ylabel('activations', fontsize=7)
    ax[1].set_xlabel('time from stimulus onset (s)')
    ax[0].set_ylabel('signal reaching the gate')
    # panel letter: this figure is panel B of the manuscript's Figure 2, whose panel A is a
    # hand-drawn schematic; the panel letter is added here rather than after export
    fig.text(0.005, 0.94, 'B', fontsize=20, weight='bold', va='top')
    fig.savefig(fname, bbox_inches='tight')
    if png:
        fig.savefig(png, dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    R = json.load(open('results/representative.json'))
    P2 = dict(R['params'], alpha=R['alpha'])
    FIG = '../manuscript/latex/figures/'
    fig_trajectories(2, fname='plots/fig_traj_m2.pdf', png=FIG + 'fig01.png')
    fig_trajectories(1, fname='plots/fig_traj_m1.pdf', png=FIG + 'fig02.png')
    fig_single(P2, png=FIG + 'fig04.png')
    print(f"figures written at the representative draw {R['draw']} "
          f"(beta = {R['beta']:.4f}, the Model 2 admissible median)")
