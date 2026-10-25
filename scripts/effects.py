"""
effects.py -- the effects compared between the models and the data.

Every effect is a fixed linear combination of the six coherence x duration condition
values of one quantity: accuracy, response time, premotor time (PMT), motor time (MT), or
the mean level, slope or curvature of the PMA rate in one analysis window. It is computed
in the same way for each participant (observed_per_subject) and for each simulated
parameter vector (draw_effects), and its Monte-Carlo standard error is propagated from the
per-condition errors the sweep records (effects_se). Used by classify.py, leak_summary.py, fig_effects.py
and table_comparison.py.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.targets import build

# Contrast weights, one per effect, applied to the six condition values in the order of
# targets.CELLS: 0-2 = 0.3 s at 0/12/50% coherence, 3-5 = 2.0 s at 0/12/50%. The name of
# each effect says which difference is taken.
DUR = np.array([1, 1, 1, -1, -1, -1]) / 3.0          # short minus long
COH = np.array([-.5, 0, .5, -.5, 0, .5])             # 50% minus 0% coherence
W = {
    # --- behavior
    'acc_coh':   COH,
    'acc_dur':  -DUR,                                 # long minus short
    'acc_int':   np.array([.5, -1, .5, -.5, 1, -.5]),  # duration benefit at 12% vs 0%, 50%
    'rt_coh':   -COH,                                 # 0% minus 50%
    'rt_dur':    DUR,                                 # short minus long
    'rt_int':    np.array([1, 0, -1, -1, 0, 1]),      # coherence effect, short minus long
    # --- stimulus-onset window
    'stim_coh':  COH,
    'stim_dur':  DUR,
    # --- delay window
    'delay_coh': COH,
    'delay_dur': -DUR,                                # long minus short
    'delay_int': np.array([1, 0, -1, -1, 0, 1]),      # coherence effect, long minus short
    'trend':     np.ones(6) / 6.0,
    'trend_dur': -DUR,                                # long minus short
    'curv':      np.ones(6) / 6.0,
    'stim_trend_coh': COH,                            # 50% minus 0%
    'stim_trend_dur': DUR,                            # short minus long
    # --- PMT and MT
    'pmet_coh':    -COH,                              # 0% minus 50%
    'pmet_dur':     DUR,                              # short minus long
    'pmet_int':     np.array([1, 0, -1, -1, 0, 1]),   # coherence effect, short minus long
    'met_dur':      DUR,                              # short minus long
}
# Which per-cell quantity each contrast is taken on, and the scale it is reported in.
SRC = {
    'acc_coh': ('acc', 1.0), 'acc_dur': ('acc', 1.0), 'acc_int': ('acc', 1.0),
    'rt_coh': ('rt', 1000.0), 'rt_dur': ('rt', 1000.0), 'rt_int': ('rt', 1000.0),
    'stim_coh': ('srate', 1.0), 'stim_dur': ('srate', 1.0),
    'delay_coh': ('rate', 1.0), 'delay_dur': ('rate', 1.0), 'delay_int': ('rate', 1.0),
    'trend': ('trend', 1.0), 'trend_dur': ('trend', 1.0), 'curv': ('curv', 1.0),
    'pmet_coh': ('pmet', 1000.0), 'pmet_dur': ('pmet', 1000.0),
    'pmet_int': ('pmet', 1000.0), 'met_dur': ('met', 1000.0),
    'stim_trend_coh': ('stim_trend', 1.0), 'stim_trend_dur': ('stim_trend', 1.0),
}
# Panel groups, in the order the Results take them. Each entry is (window, [keys]); the
# figure draws one labelled block per window, three panels to a row.
GROUPS = [
    ('Accuracy, RT, PMT and MT', ['acc_coh', 'acc_dur', 'acc_int', 'rt_coh', 'rt_dur',
                                  'rt_int', 'pmet_coh', 'pmet_dur', 'pmet_int', 'met_dur']),
    ('Stimulus-onset window', ['stim_coh', 'stim_dur', 'stim_trend_coh', 'stim_trend_dur']),
    ('Delay window', ['delay_coh', 'delay_dur', 'delay_int', 'trend', 'trend_dur', 'curv']),
]
# Effects whose direction the derivations fix at every parameter setting (Results,
# Theoretical predictions): no sign summary is reported for them (Materials and Methods,
# Model-data comparison).
SIGN_FIXED = {'acc_coh', 'acc_dur', 'rt_coh', 'rt_dur', 'pmet_coh', 'pmet_dur',
              'stim_coh', 'delay_coh', 'delay_dur', 'trend'}
KEYS = [k for _, ks in GROUPS for k in ks]

NAMES = dict(
    acc_coh='accuracy,\n50% minus 0% coherence',
    acc_dur='accuracy,\nlong minus short',
    acc_int='duration benefit on accuracy,\nat 12% vs 0% and 50%',
    rt_coh='RT,\n0% minus 50% coherence (ms)',
    rt_dur='RT, short minus long (ms)',
    rt_int='coherence effect on RT,\nshort minus long (ms)',
    pmet_coh='PMT,\n0% minus 50% coherence (ms)',
    pmet_dur='PMT, short minus long (ms)',
    pmet_int='coherence effect on PMT,\nshort minus long (ms)',
    met_dur='MT, short minus long (ms)',
    stim_coh='mean PMA rate,\n50% minus 0% coherence (/s)',
    stim_dur='mean PMA rate,\nshort minus long (/s)',
    stim_trend_coh='PMA-rate slope,\n50% minus 0% coherence (/s$^2$)',
    stim_trend_dur='PMA-rate slope,\nshort minus long (/s$^2$)',
    delay_coh='mean PMA rate,\n50% minus 0% coherence (/s)',
    delay_dur='mean PMA rate,\nlong minus short (/s)',
    delay_int='coherence effect on mean PMA rate,\nlong minus short (/s)',
    trend='PMA-rate slope (/s$^2$)',
    trend_dur='PMA-rate slope,\nlong minus short (/s$^2$)',
    curv='PMA-rate curvature (/s$^3$)')


def effects(cells):
    """`cells` maps a per-cell quantity name to its six values."""
    return {k: float(np.asarray(W[k], float) @ np.asarray(cells[SRC[k][0]], float)
                     * SRC[k][1]) for k in KEYS}


def effects_se(m):
    """Monte-Carlo standard error of each effect, propagated from the per-condition errors
    the sweep records: binomial for accuracy, the between-trial error of a mean for RT,
    PMT and MT, and the Poisson error of the event counts for the PMA rates, slopes and
    curvature. Count-based errors are floored at one event, so that a condition with no
    event is not credited with zero uncertainty."""
    g = lambda k: np.asarray(m['raw'][k], float)
    acc, n = g('acc'), g('rt_n')

    def bse(p, ntr):
        p = np.maximum(p, 1.0 / ntr)
        return np.sqrt(np.clip(p * (1 - p), 0, None) / ntr)

    def pse(rate, expo):
        return np.sqrt(np.maximum(rate * expo, 1.0)) / np.maximum(expo, 1)

    se = dict(acc=bse(acc, np.maximum(n, 1)),
              rt=g('rt_se'), pmet=g('pmet_se'), met=g('met_se'),
              trend=g('trend_se'), curv=g('curv_se'),
              stim_trend=g('stim_trend_se'),
              srate=pse(g('srate'), n * 2.0), rate=pse(g('rate'), n * 2.0))
    out = {}
    for k in KEYS:
        src, scale = SRC[k]
        w = np.asarray(W[k], float)
        out[k] = float(np.sqrt((w ** 2) @ (np.asarray(se[src], float) ** 2))) * scale
    return out


def observed_per_subject():
    """One dict of effects per participant, computed from the data exactly as from a
    simulated draw."""
    E = build()
    nsub = E['rt'].shape[0]
    obs_cells = [dict(acc=E['acc'][s], rt=E['rt'][s],
                      srate=E['stim'][s].sum(axis=-1) / 2.0,
                      rate=E['dly'][s].sum(axis=-1) / 2.0,
                      trend=E['trend'][s], curv=E['quad'][s],
                      stim_trend=E['stim_trend'][s],
                      pmet=E['pmet'][s], met=E['met'][s]) for s in range(nsub)]
    return [effects(c) for c in obs_cells]


def draw_effects(m):
    """The effects of one simulated (draw, model) record."""
    return effects({k: np.asarray(m['raw'][src], float) for k, src in
                    (('acc', 'acc'), ('rt', 'rt_mean'), ('srate', 'srate'),
                     ('rate', 'rate'), ('trend', 'trend'), ('curv', 'curv'),
                     ('stim_trend', 'stim_trend'),
                     ('pmet', 'pmet_mean'), ('met', 'met_mean'))})
