"""
figures_space.py -- the two figures that come from the parameter-space evaluation rather
than from a single simulated trial.

  fig_frontier.pdf       Model 1's ramp-rate frontier, with the sampled draws on it. This
                         is the figure that carries the delay-period argument: the two
                         quantities that have to be reproduced together are on the two
                         axes, the analytic frontier is the boundary of what Model 1 can
                         do, and the observation sits far outside it.
  fig_delay_profile.pdf  The delay-period rate curve itself, observed against what each
                         architecture produces over its admissible space.

Inputs: results/sweep_screened.jsonl, results/frontier.json, results/ramp_check.jsonl,
pma_counts.csv.  Usage: python3 figures_space.py
"""
import json, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
C_M1, C_M2, C_OBS = '#c1272d', '#3182bd', '#222222'

# Observed delay-period statistics, from the 132 participant x condition cells.
# OBS_RATE is the MEAN rate over the 1-3 s window -- the quantity the simulation reports,
# and NOT the end-of-delay peak (1.446 s^-1), which is a different number.
OBS_RATE, OBS_RMIN, OBS_RMAX = 0.709, 0.085, 2.785
OBS_SLOPE, OBS_SLO, OBS_SHI = 0.817, 0.659, 1.029

recs = [json.loads(l) for l in open('results/sweep_screened.jsonl')]
FR = json.load(open('results/frontier.json'))


def per_model(tag):
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    rate = np.array([np.nanmean(np.asarray(m['raw']['rate'], float)) for m in adm])
    slope = np.array([np.nanmean(np.asarray(m['raw']['slope'], float)) for m in adm])
    ok = np.isfinite(rate) & np.isfinite(slope) & (rate > 0)
    prof = np.array([np.asarray(m['raw']['prof'], float).sum(axis=0) for m in adm])
    expo = np.array([np.asarray(m['raw']['prof_expo'], float).sum() for m in adm])
    return rate[ok], slope[ok], prof[ok], expo[ok]


def observed_profile():
    cnt = np.zeros(20); expo = np.zeros(20)
    for r in csv.DictReader(open('pma_counts.csv')):
        if r['period'] != 'delay':
            continue
        b = int(r['bin']); cnt[b] += float(r['count']); expo[b] += float(r['ntrial'])
    return cnt, expo


# ------------------------------------------------------------------ 1. the frontier
def fig_frontier(fname='plots/fig_frontier.pdf'):
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    env = sorted(FR['envelope'], key=lambda r: r['rate'])
    x = np.array([e['rate'] for e in env]); y = np.array([e['slope_at_least'] for e in env])
    # shade what Model 1 CANNOT reach, not what it can: the eye reads a shaded region as
    # excluded, and the excluded region is the one above the frontier
    ax.fill_between(x, y, 4.0, color=C_M1, alpha=.06, lw=0, zorder=0)
    ax.plot(x, y, color=C_M1, lw=1.6, zorder=3,
            label='Model 1: largest attainable slope (exact)')
    ax.text(0.011, 2.15, 'beyond Model 1', color=C_M1, fontsize=7.5, style='italic')
    for tag, col, lab in [('model1', C_M1, 'Model 1, admissible draws'),
                          ('model2', C_M2, 'Model 2, admissible draws')]:
        r, s, _, _ = per_model(tag)
        ax.scatter(r, s, s=9, facecolor='none', edgecolor=col, lw=.6, alpha=.75,
                   zorder=2, label=lab)
    rc = [json.loads(l) for l in open('results/ramp_check.jsonl')]
    by = {}
    for d in rc:
        by.setdefault(d['draw'], []).append(d)
    R = np.array([np.mean([x['rate'] for x in v]) for v in by.values()])
    S = np.array([np.mean([x['slope'] for x in v]) for v in by.values()])
    ax.scatter(R, S, s=26, marker='D', color=C_M1, zorder=4,
               label='Model 1, re-simulated at 40 000 trials')
    ax.axhspan(OBS_SLO, OBS_SHI, color='0.86', zorder=1)
    ax.errorbar([OBS_RATE], [OBS_SLOPE],
                xerr=[[OBS_RATE - OBS_RMIN], [OBS_RMAX - OBS_RATE]], fmt='none',
                ecolor='0.55', elinewidth=.8, capsize=2, zorder=5)
    ax.plot([OBS_RATE], [OBS_SLOPE], marker='*', ms=13, color=C_OBS, zorder=6,
            label='observed')
    ax.set_xscale('log')
    ax.set_xlim(8e-3, 4.0); ax.set_ylim(-1.7, 2.6)
    ax.set_xlabel('delay-period partial-activation rate (s$^{-1}$)')
    ax.set_ylabel('delay-period log-rate slope (s$^{-1}$)')
    ax.axhline(0, color='0.75', lw=.7, ls=':')
    ax.legend(loc='lower left', frameon=False, handletextpad=.5)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname)


# ------------------------------------------------------------ 2. the delay profile
def fig_delay_profile(fname='plots/fig_delay_profile.pdf', lo=0.4, hi=1.2):
    """The delay profile, at draws that produce activations at about the observed rate.

    Conditioning on the rate is what makes the picture readable and is also the honest
    question: over the whole admissible set the levels differ by two orders of magnitude,
    so a percentile envelope taken bin by bin mixes profiles that are not comparable and
    shows nothing about shape. Restricting to draws whose mean rate brackets the observed
    0.709 s-1 asks the question the delay argument actually turns on -- when the
    architecture produces partial activations as often as the participants do, does its
    rate curve rise as theirs does?
    """
    t = np.arange(20) * 0.1 + 0.05
    cnt, expo = observed_profile()
    obs = cnt / (expo * 0.1)
    se = obs / np.sqrt(np.maximum(cnt, 1))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True,
                           gridspec_kw=dict(wspace=0.08))
    for a, (tag, col, name) in zip(ax, [('model1', C_M1, 'Model 1'),
                                        ('model2', C_M2, 'Model 2')]):
        r, _, prof, pex = per_model(tag)
        sel = (r >= lo) & (r <= hi)
        rate = prof[sel] / (pex[sel][:, None] * 0.1)
        band = np.percentile(rate, [5, 95], axis=0)
        a.fill_between(t, band[0], band[1], color=col, alpha=.15, lw=0,
                       label=f'{name}, 5th–95th percentile ({sel.sum()} draws)')
        a.plot(t, np.median(rate, axis=0), color=col, lw=1.6, label=f'{name}, median')
        a.errorbar(t, obs, yerr=se, fmt='o', ms=2.6, color=C_OBS, lw=.8, label='observed')
        a.set_xlabel('time into the delay analysis window (s)')
        a.legend(loc='upper left', frameon=False)
        a.set_ylim(0, 2.0)
    ax[0].set_ylabel('partial-activation rate (s$^{-1}$)')
    ax[0].set_title(f'draws with a mean rate in [{lo}, {hi}] s$^{{-1}}$ '
                    f'(observed {OBS_RATE})', fontsize=8, loc='left')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname)


if __name__ == '__main__':
    os.makedirs('plots', exist_ok=True)
    fig_frontier()
    fig_delay_profile()
