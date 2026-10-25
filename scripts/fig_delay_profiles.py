"""
fig_delay_profiles.py -- the shape of the delay-window PMA rate (manuscript Figure 5):
one line per admissible parameter vector of each model, each divided by its own mean over
the delay analysis window, against the observed profile normalised the same way.

Reads results/sweep_screened.jsonl and the archive (DELAY_DATA). Writes
plots/fig_delay_profiles.pdf and the manuscript's copy, ../manuscript/latex/figures/fig07.png.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.targets import build

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
C_M1, C_M2, C_OBS = '#c1272d', '#3182bd', '#222222'
PNG = '../manuscript/latex/figures/fig07.png'
recs = [json.loads(l) for l in open('results/sweep_screened.jsonl')]


def per_model(tag, stat='trend'):
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    rate = np.array([np.nanmean(np.asarray(m['raw']['rate'], float)) for m in adm])
    y = np.array([np.nanmean(np.asarray(m['raw'][stat], float)) for m in adm])
    ok = np.isfinite(rate) & np.isfinite(y) & (rate > 0)
    prof = np.array([np.asarray(m['raw']['prof'], float).sum(axis=0) for m in adm])
    expo = np.array([np.asarray(m['raw']['prof_expo'], float).sum() for m in adm])
    return rate[ok], y[ok], prof[ok], expo[ok]



def fig_delay_profiles(fname=None):
    """One line per admissible parameter vector: its delay-window PMA rate divided by its
    own mean over the window, against the observed profile normalised the same way. A
    vector enters if its profile rests on at least 200 events."""
    fname = fname or 'plots/fig_delay_profiles.pdf'
    t = np.arange(20) * 0.1 + 0.05
    obs = build()['dly'].mean(axis=(0, 1))          # per-trial counts per bin, averaged
    obs_n = obs / obs.mean()
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True,
                           gridspec_kw=dict(wspace=0.08))
    for a, (tag, col, name) in zip(ax, [('model1', C_M1, 'Model 1'),
                                        ('model2', C_M2, 'Model 2')]):
        _, _, prof, pex = per_model(tag)
        rate = prof / (pex[:, None] * 0.1)
        m = rate.mean(axis=1)
        keep = (m > 0) & (prof.sum(axis=1) >= 200)
        norm = rate[keep] / m[keep][:, None]
        for y in norm:
            a.plot(t, y, color=col, lw=.4, alpha=.30, zorder=1)
        a.plot(t, np.median(norm, axis=0), color=col, lw=2.0, zorder=3)
        a.plot(t, obs_n, 'o-', ms=3, color=C_OBS, lw=1.2, zorder=4)
        a.set_xlabel('time into the delay analysis window (s)')
        a.set_title(f'{name}, {keep.sum()} admissible parameter vectors', fontsize=8)
        a.set_ylim(0, 3.2)
        a.axhline(1.0, color='0.75', lw=.7, ls=':')
    ax[0].set_ylabel('PMA rate,\nrelative to its own window mean')
    ax[1].legend([plt.Line2D([], [], color=C_M2, lw=.4),
                  plt.Line2D([], [], color=C_M2, lw=2.0),
                  plt.Line2D([], [], color=C_OBS, lw=1.2, marker='o', ms=3)],
                 ['one admissible parameter vector', 'median over vectors', 'observed'],
                 fontsize=6.5, frameon=False, loc='upper left')
    fig.savefig(fname, bbox_inches='tight')
    fig.savefig(PNG, dpi=200, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname, 'and', PNG)


if __name__ == '__main__':
    os.makedirs('plots', exist_ok=True)
    fig_delay_profiles()
