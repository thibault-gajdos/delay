"""
fig_effects.py -- the model-data comparison, effect by effect (manuscript Figure 8).

For every effect and each model: the range of the effect over the admissible parameter
vectors (minimum to maximum), against the observed group mean and its 95% confidence
interval, with the number of vectors that reproduce the observed magnitude, do not
reproduce it, or are numerically uncertain (classify.py). Panels are grouped as the
Results take them (effects.GROUPS).

Reads results/classification.json. Writes plots/fig_effects.pdf and the manuscript's copy,
../manuscript/latex/figures/fig08.png.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from effects import GROUPS, NAMES

C = json.load(open('results/classification.json'))
COL = {'model1': '#c1272d', 'model2': '#3182bd'}
LAB = {'model1': 'Model 1', 'model2': 'Model 2'}
NCOL = 4
plt.rcParams.update({'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.spines.left': False})

rows = []                 # ('head', group name) or ('panels', keys), top to bottom
for name, keys in GROUPS:
    rows.append(('head', name))
    for i in range(0, len(keys), NCOL):
        rows.append(('panels', keys[i:i + NCOL]))
heights = [0.32 if kind == 'head' else 1.0 for kind, _ in rows]

fig = plt.figure(figsize=(7.2, 1.3 * sum(heights) + 0.4))
gs = fig.add_gridspec(len(rows), NCOL, height_ratios=heights, hspace=0.9, wspace=0.35,
                      left=0.03, right=0.97, top=0.99, bottom=0.04)
letter = 0
for r, (kind, what) in enumerate(rows):
    if kind == 'head':
        ax = fig.add_subplot(gs[r, :]); ax.axis('off')
        ax.text(0, 0, what, fontsize=8.5, weight='bold', va='bottom', transform=ax.transAxes)
        continue
    for c, k in enumerate(what):
        ax = fig.add_subplot(gs[r, c])
        o = C['observed'][k]
        ax.axvspan(o['lo'], o['hi'], color='0.85', zorder=0)
        ax.axvline(o['mean'], color='k', lw=0.8, zorder=1)
        ax.axvline(0, color='0.6', lw=0.5, ls=':', zorder=1)
        for y, tag in ((1, 'model1'), (0, 'model2')):
            m = C['models'][tag][k]
            ax.plot([m['min'], m['max']], [y, y], color=COL[tag], lw=4,
                    solid_capstyle='butt', zorder=2)
            ax.text(1.0, y, f"{m['reproduces']}/{m['does_not']}/{m['uncertain']}",
                    transform=ax.get_yaxis_transform(), fontsize=5.5, va='center',
                    ha='left', color=COL[tag])
        ax.set_ylim(-0.7, 1.7); ax.set_yticks([])
        ax.tick_params(axis='x', labelsize=6, length=2)
        ax.set_title(NAMES[k], fontsize=6.5, loc='left')
        ax.text(-0.04, 1.32, chr(ord('A') + letter), transform=ax.transAxes,
                fontsize=8, weight='bold', va='top', ha='right')
        letter += 1
handles = [plt.Line2D([], [], color=COL[t], lw=4, label=LAB[t]) for t in ('model1', 'model2')]
handles += [plt.Line2D([], [], color='k', lw=0.8, label='observed mean'),
            plt.Rectangle((0, 0), 1, 1, color='0.85', label='observed 95% CI')]
fig.legend(handles=handles, loc='lower right', fontsize=6.5, frameon=False, ncol=2)
os.makedirs('plots', exist_ok=True)
fig.savefig('plots/fig_effects.pdf', bbox_inches='tight')
fig.savefig('../manuscript/latex/figures/fig08.png', dpi=200, bbox_inches='tight')
print('wrote plots/fig_effects.pdf and ../manuscript/latex/figures/fig08.png')
