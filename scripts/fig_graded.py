"""
fig_graded.py -- manuscript Figure S2: graded maintenance at the representative parameter
vector. Mean signal reaching the gate during the delay, in units of the internal threshold
theta, for each coherence, after the 0.3-s and the 2-s stimulus, under Model 1 and Model 2.
Writes ../manuscript/latex/figures/fig11.png and plots/fig_graded.pdf.
Adapted from J. Achard's fig_graded_supp.py.

The two models are simulated with the same random seeds, so that they differ only in beta.
Trials are simulated in chunks to bound memory. Run representative.py first.
Usage: python3 fig_graded.py
"""
import numpy as np, json, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS

NTRIALS, CHUNK = 10000, 2000
TG = np.linspace(0.0, 3.0, 350)            # time from physical stimulus offset (s)
DURS, COHS = [0.3, 2.0], [0.0, 0.12, 0.50]
COL = {0.3: ['#E5B1EF', '#B23FC7', '#76248C'],   # purple, light to dark = 0, 12, 50%
       2.0: ['#A1D9B0', '#258040', '#104F29']}   # green
plt.rcParams.update({'font.size': 8, 'lines.linewidth': 1.4})


def mean_curve(P, model, c, D, k):
    """Mean over valid-response trials of the signal reaching the gate, on the grid TG."""
    total, count = np.zeros_like(TG), np.zeros_like(TG)
    for b in range(NTRIALS // CHUNK):
        s = simulate(P, c, D, CHUNK, model, seed=5000 + 100 * k + b)
        v = np.flatnonzero(measures(s)['valid'])
        Y = s['Z'][v] + P['g']
        t = np.arange(Y.shape[1]) * DT - D
        M = np.array([np.interp(TG, t, y, left=np.nan, right=np.nan) for y in Y])
        ok = np.isfinite(M)
        total += np.where(ok, M, 0.0).sum(0)
        count += ok.sum(0)
    return total / np.maximum(count, 1)


if __name__ == '__main__':
    R = json.load(open('results/representative.json'))
    theta = float(R['params']['theta'])
    P = {1: dict(R['params'], alpha=0.0), 2: dict(R['params'], alpha=float(R['alpha']))}
    curves = {}
    for model in (1, 2):
        for k, (c, D) in enumerate(CELLS):
            curves[(model, D, c)] = mean_curve(P[model], model, c, D, k) / theta
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey='row')
    for row, D in enumerate(DURS):
        for col, model in enumerate((1, 2)):
            ax = axes[row, col]
            for j, c in enumerate(COHS):
                ax.plot(TG, curves[(model, D, c)], color=COL[D][j], lw=2,
                        label=f'{round(100 * c)}%')
            ax.set_xlim(0, 3)
            ax.set_title(f'Model {model}, {D:g}-s stimulus', fontsize=12)
            ax.grid(axis='y', color='#D0D0D0', lw=0.6, alpha=0.7)
            if row == 1:
                ax.set_xlabel('Time from stimulus offset (s)', fontsize=10)
            if col == 0:
                ax.set_ylabel(r'signal reaching the gate / $\theta$', fontsize=10)
                ax.legend(title='Coherence', loc=('lower right' if row == 0 else 'upper right'), frameon=False,
                          fontsize=10, title_fontsize=10)
            ax.text(-0.12, 1.04, 'ABCD'[2 * row + col], transform=ax.transAxes,
                    fontsize=12, weight='bold', va='bottom')
    fig.tight_layout()
    fig.savefig('plots/fig_graded.pdf', bbox_inches='tight')
    fig.savefig('../manuscript/latex/figures/fig11.png', dpi=200, bbox_inches='tight')
    print('wrote ../manuscript/latex/figures/fig11.png')
