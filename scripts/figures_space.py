"""
figures_space.py -- the two figures that come from the parameter-space evaluation rather
than from a single simulated trial.

  fig_rate_trend.pdf     The figure that carries the delay-period argument: the two
                         quantities that have to be reproduced together -- the mean rate
                         over the delay window and its linear trend -- on the two axes,
                         with the admissible draws of both architectures and the
                         observation. Under the stationary start Model 1's trend is zero
                         analytically, so the panel needs no computed frontier: the
                         boundary of what Model 1 can do is the horizontal axis.
  fig_delay_profiles.pdf The delay-period rate curve itself, observed against what each
                         architecture produces over its admissible space.

Inputs: results/sweep_screened.jsonl, results/effect_panel.json, pma_counts.csv.
Usage: python3 figures_space.py [DUMP.jsonl]
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
# OBS_TREND is the delay-period linear contrast of the RAW rate (Methods: TREND), in s^-2,
# read from the archive's own bin export so that it is the same number the Results report.
# Being linear in the rate it has no pooled-versus-per-cell ambiguity: the mean of the
# 132 participant x condition contrasts and the contrast of the pooled profile agree at
# +0.5604. OBS_TLO/OBS_THI span the six CELL MEANS -- a spread across conditions, not an
# uncertainty -- and OBS_TCI is the between-participant 95% confidence interval, which is
# the interval to compare a model against.
OBS_TREND, OBS_TLO, OBS_THI = 0.5604, 0.4685, 0.6342
OBS_TCI = (0.3753, 0.7454)

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
recs = [json.loads(l) for l in open(DUMP)]


def per_model(tag, stat='trend'):
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    rate = np.array([np.nanmean(np.asarray(m['raw']['rate'], float)) for m in adm])
    y = np.array([np.nanmean(np.asarray(m['raw'][stat], float)) for m in adm])
    ok = np.isfinite(rate) & np.isfinite(y) & (rate > 0)
    prof = np.array([np.asarray(m['raw']['prof'], float).sum(axis=0) for m in adm])
    expo = np.array([np.asarray(m['raw']['prof_expo'], float).sum() for m in adm])
    return rate[ok], y[ok], prof[ok], expo[ok]


def observed_profile():
    cnt = np.zeros(20); expo = np.zeros(20)
    for r in csv.DictReader(open('pma_counts.csv')):
        if r['period'] != 'delay':
            continue
        b = int(r['bin']); cnt[b] += float(r['count']); expo[b] += float(r['ntrial'])
    return cnt, expo


# ------------------------------------------- 1. rate against trend, both architectures
def fig_rate_trend(fname=None):
    """The delay-period argument in one panel: the window-mean rate on the x axis, its
    linear trend on the y axis, one point per admissible parameter set.

    The two quantities have to be reproduced TOGETHER, which is the whole reason the trend
    replaced the log-rate slope: a draw can earn a steep log-rate slope by ramping from
    almost nothing, and the raw-rate trend cannot be earned that way because it is the
    product of the level and the relative ramp. Model 1's trend is zero analytically under
    the stationary start, so its cloud has to sit on the dotted line whatever its rate --
    no envelope is computed here, and none is needed.
    """
    fname = fname or 'plots/fig_rate_trend.pdf'
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    for tag, col, lab in [('model1', C_M1, 'Model 1, admissible draws'),
                          ('model2', C_M2, 'Model 2, admissible draws')]:
        r, y, _, _ = per_model(tag)
        ax.scatter(r, y, s=11, facecolor='none', edgecolor=col, lw=.6, alpha=.75,
                   zorder=2, label=lab)
    ax.axhline(0, color=C_M1, lw=1.4, ls=(0, (4, 2)), zorder=3)
    ax.set_ylim(-0.26, 1.26)
    # the caption for Model 1's line goes BELOW it, in the band the draws leave empty:
    # placed above it, it sat on top of Model 1's own cloud.
    # The simulated points sit slightly BELOW the analytic line, and the offset is a
    # measurement bias rather than a prediction: the response-related burst is excluded
    # from the event count wherever it starts, and in the model it starts inside the last
    # 200 ms of the delay window on up to 2% of trials, against 0.02% in the data. It
    # biases both architectures' trends downward by about 0.01 s^-2, so Model 2's numbers
    # are conservative and Model 1's departure from zero is entirely accounted for.
    ax.text(0.0095, -0.245, 'Model 1: zero trend at every admissible setting (largest '
                            'attainable over the\nsampled ranges, $3\\times10^{-11}$ '
                            's$^{-2}$). The scatter of the simulated points is\nMonte-Carlo noise and their small negative offset is the excluded\nresponse burst; neither is a prediction.',
            color=C_M1, fontsize=7.0, va='bottom', linespacing=1.5)
    ax.axhspan(*OBS_TCI, color='0.86', zorder=0)
    ax.errorbar([OBS_RATE], [OBS_TREND],
                xerr=[[OBS_RATE - OBS_RMIN], [OBS_RMAX - OBS_RATE]], fmt='none',
                ecolor='0.55', elinewidth=.8, capsize=2, zorder=5)
    ax.plot([OBS_RATE], [OBS_TREND], marker='*', ms=13, color=C_OBS, zorder=6,
            label='observed (95% CI shaded)')
    ax.set_xscale('log')
    ax.set_xlim(8e-3, 4.0)
    ax.set_xlabel('delay-period partial-activation rate (s$^{-1}$)')
    ax.set_ylabel('delay-period rate trend (s$^{-2}$)')
    ax.legend(loc='upper left', frameon=False, handletextpad=.5)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname)


# ------------------------------------------------------------ 2. the delay profile
def fig_delay_profiles(fname=None):
    """The shape of the delay-period rate curve, one line per admissible parameter set.

    Each thin line is the delay-period profile of one admissible draw, divided by its own
    mean over the window so that draws whose rates differ by two orders of magnitude can be
    compared on shape; the observed profile is normalised the same way. Showing the family
    rather than a median or a percentile band is the point: the question is whether the
    architecture CAN produce the observed rise, not what it does at a typical setting, and a
    band taken bin by bin mixes profiles that are not comparable.
    """
    fname = fname or 'plots/fig_delay_profiles.pdf'
    t = np.arange(20) * 0.1 + 0.05
    cnt, expo = observed_profile()
    obs = cnt / (expo * 0.1)
    obs_n = obs / obs.mean()
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True,
                           gridspec_kw=dict(wspace=0.08))
    for a, (tag, col, name) in zip(ax, [('model1', C_M1, 'Model 1'),
                                        ('model2', C_M2, 'Model 2')]):
        _, _, prof, pex = per_model(tag)
        rate = prof / (pex[:, None] * 0.1)
        m = rate.mean(axis=1)
        # a draw enters only if its profile rests on at least 200 events over the window:
        # below that the normalised curve is Poisson noise rather than shape. This is a
        # precision criterion on the estimate, not a selection on what the estimate says.
        keep = (m > 0) & (prof.sum(axis=1) >= 200)
        norm = rate[keep] / m[keep][:, None]
        for y in norm:
            a.plot(t, y, color=col, lw=.4, alpha=.30, zorder=1)
        a.plot(t, np.median(norm, axis=0), color=col, lw=2.0, zorder=3)
        a.plot(t, obs_n, 'o-', ms=3, color=C_OBS, lw=1.2, zorder=4)
        a.set_xlabel('time into the delay analysis window (s)')
        a.set_title(f'{name}, {keep.sum()} admissible parameter sets', fontsize=8)
        a.set_ylim(0, 3.2)
        a.axhline(1.0, color='0.75', lw=.7, ls=':')
    ax[0].set_ylabel('partial-activation rate,\nrelative to its own window mean')
    ax[1].legend([plt.Line2D([], [], color=C_M2, lw=.4),
                  plt.Line2D([], [], color=C_M2, lw=2.0),
                  plt.Line2D([], [], color=C_OBS, lw=1.2, marker='o', ms=3)],
                 ['one admissible parameter set', 'median over parameter sets', 'observed'],
                 fontsize=6.5, frameon=False, loc='upper left')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname)


# ---------------------------------------------------- 3. the effect-size panel
def fig_effect_sizes(fname=None):
    """What each model can produce, on the effects that carry an argument, against
    the observation. Bars span the 5th-95th percentile over admissible draws, the tick is
    the median; the vertical line is the observed value with its 95% CI shaded. A bar that
    does not reach the line is a quantitative failure of that model on that effect.
    Panels are lettered A-I in the order of effect_panel.NAMES."""
    fname = fname or 'plots/fig_effect_sizes.pdf'
    P = json.load(open('results/effect_panel.json'))
    keys = list(P['names'])
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(8.4, 2.0 * nrow))
    axes = np.atleast_2d(axes).ravel()
    for a in axes[len(keys):]:
        a.axis('off')
    for a, k in zip(axes, keys):
        for y, (tag, col) in enumerate([('model2', C_M2), ('model1', C_M1)]):
            r = P['models'][tag][k]
            a.plot([r['p5'], r['p95']], [y, y], color=col, lw=4, solid_capstyle='butt',
                   alpha=.85, zorder=2)
            if not r['covers_observed']:
                a.plot([r['p95'] if r['p95'] < P['observed'][k] else r['p5']], [y],
                       marker='x', ms=5, color=col, mew=1.4, zorder=4)
        o, se = P['observed'][k], P['observed_se'][k]
        a.axvspan(o - 1.96 * se, o + 1.96 * se, color='0.86', zorder=0)
        a.axvline(o, color=C_OBS, lw=1.1, zorder=1)
        a.set_yticks([0, 1]); a.set_yticklabels([])
        a.set_ylim(-0.6, 1.6)
        a.set_title(P['names'][k], fontsize=7.0, pad=4)
        # panel letter, in the order of NAMES (effect_panel.py); the Results cite these
        a.text(-0.30, 1.22, chr(65 + list(keys).index(k)), transform=a.transAxes,
               fontsize=9.5, fontweight='bold', va='top', ha='left')
        a.tick_params(axis='x', labelsize=6.5)
        # the DIRECTION, printed per model as two percentages: the share of admissible
        # sets whose effect is, at two Monte-Carlo standard errors, significantly of the
        # observed sign and the share significantly against it. The remainder is not
        # distinguishable from zero and carries no directional information -- separating
        # it out matters, because on the premature-response effect more than half of the
        # sets produce too few premature responses to test the prediction at all, and on
        # Model 1's delay trend the sets that are not negative are simply flat.
        for y, tag in enumerate(['model2', 'model1']):
            r = P['models'][tag][k]
            col = C_M2 if tag == 'model2' else C_M1
            a.text(1.03, y, f"{r['frac_right']*100:.0f}", transform=a.get_yaxis_transform(),
                   va='center', fontsize=6.5, color=col)
            a.text(1.16, y, '/', transform=a.get_yaxis_transform(), va='center',
                   fontsize=6.5, color='0.6')
            a.text(1.22, y, f"{r['frac_wrong']*100:.0f}", transform=a.get_yaxis_transform(),
                   va='center', fontsize=6.5, color=col,
                   fontweight=('bold' if r['frac_wrong'] > 0.25 else 'normal'))
    for a in axes[len(keys):]:
        a.set_visible(False)
    for r in range(nrow):
        axes[r * ncol].set_yticklabels(['Model 2', 'Model 1'], fontsize=7.5)
    axes[0].plot([], [], color=C_OBS, lw=1.1, label='observed (95% CI shaded)')
    fig.subplots_adjust(wspace=0.62, hspace=0.95, bottom=0.16)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower left',
               bbox_to_anchor=(0.06, 0.045), frameon=False, fontsize=7)
    fig.text(0.06, 0.005, 'the two figures at the right of each bar are the percentages of '
             'admissible parameter sets whose effect is, at two Monte-Carlo standard errors, '
             'significantly of the observed sign / significantly against it;\nthe remainder '
             'is not distinguishable from zero. A cross marks a model whose 5th-95th '
             'percentile interval does not reach the observed value.',
             fontsize=6.6, va='bottom', color='0.35')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)
    print('wrote', fname)


if __name__ == '__main__':
    os.makedirs('plots', exist_ok=True)
    fig_effect_sizes()
    fig_rate_trend()
    fig_delay_profiles()
