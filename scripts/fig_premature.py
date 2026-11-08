"""fig_premature.py -- proportion of responses made before the response signal (RS), by
condition (manuscript Figure S1), drawn from the archive's per-participant behavioural files
(data/behavioral/s1.npy ... s22.npy). Writes ../manuscript/latex/figures/fig10.png.
Needs DELAY_DATA. Adapted from J. Achard's fig_premature_response_RS.py.

A response before the RS is a trial whose response moment (column 10) is 0 ("too early").
The proportion is taken over the 160 trials of each coherence x duration condition. Error
bars are +-1 within-subjects standard error (Cousineau-Morey normalisation), as in
fig_behavior.py.
"""
import os
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

D = os.path.join(os.environ.get('DELAY_DATA',
                 os.path.expanduser('~/thib/projects/delay/data/86cbw-osfstorage-archive')),
                 'data', 'behavioral')
NSUB = 22
STIM = {0.3: 18, 2.0: 120}            # stimulus codes in column 6
COH = [0.0, 0.12, 0.50]
COL = {0.3: ['#E5B1EF', '#B23FC7', '#76248C'],   # 0.3 s: purple, light to dark = 0, 12, 50%
       2.0: ['#A1D9B0', '#258040', '#104F29']}   # 2 s: green
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})


def proportions():
    """(participants, durations, coherences) proportion of responses before the RS."""
    P = np.zeros((NSUB, 2, 3))
    for s in range(NSUB):
        d = np.load(os.path.join(D, f's{s + 1}.npy'), allow_pickle=True)
        coh = d[:, 5].astype(float)
        stim = d[:, 6].astype(float).astype(int)
        early = d[:, 10].astype(float).astype(int) == 0
        for j, dur in enumerate(STIM):
            for k, c in enumerate(COH):
                m = (stim == STIM[dur]) & np.isclose(coh, c)
                P[s, j, k] = early[m].mean()
    return P


def morey_sem(A):
    """Within-subjects SEM: remove each participant's mean, add the grand mean, Morey correction."""
    n, M = A.shape[0], A.shape[1] * A.shape[2]
    flat = A.reshape(n, -1)
    N = (flat - flat.mean(axis=1, keepdims=True) + flat.mean()).reshape(A.shape)
    return N.std(axis=0, ddof=1) / np.sqrt(n) * np.sqrt(M / (M - 1))


P = proportions()
mean, se = P.mean(axis=0), morey_sem(P)
fig, ax = plt.subplots(figsize=(8, 6))
x = 0
for j, dur in enumerate(STIM):
    for k in range(3):
        ax.bar(x, mean[j, k], yerr=se[j, k], color=COL[dur][k], capsize=5)
        x += 1
ax.set_xticks([])
ax.set_xlabel('Conditions', fontsize=15)
ax.set_ylabel('Proportion of responses before the RS', fontsize=15)
ax.tick_params(axis='y', labelsize=14)
handles = [Patch(color=COL[dur][k]) for dur in STIM for k in (2, 1, 0)]
labels = [f'{dur:g} s, {round(100 * COH[k])}%' for dur in STIM for k in (2, 1, 0)]
ax.legend(handles, labels, title='Stimulus duration, coherence', ncol=2, frameon=True,
          loc='upper left', fontsize=10, title_fontsize=11)
fig.tight_layout()
fig.savefig('../manuscript/latex/figures/fig10.png', dpi=200, bbox_inches='tight')
print('wrote ../manuscript/latex/figures/fig10.png')
print('mean proportion per condition (rows 0.3 s, 2 s; columns 0, 12, 50%):')
print(np.round(mean, 4))
