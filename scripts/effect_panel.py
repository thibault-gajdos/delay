"""
effect_panel.py -- the effect sizes each architecture can produce, against the observed
value, for every effect on behavior or partial-activation rate that the Results report.

This is what the parameter-space evaluation is for. The predictions of the Appendix are
ordinal and hold at every admissible setting; what a derivation cannot supply is the scale.
For each effect below, this reports the distribution of that effect over a model's
admissible draws -- median and 5th-95th percentile -- next to the observed value and its
between-participant standard error, and the fraction of draws whose effect is, at two
Monte-Carlo standard errors, significantly of the observed sign, significantly against it,
or not distinguishable from zero.

The effects are grouped by ANALYSIS WINDOW, in the order the Results take them, which is
how the figure is drawn (figures_space.fig_effect_sizes). Until 2026-09-08 the panel held
only the nine effects that discriminate between the architectures, and left out the main
effects both produce; the Results, however, report those too, and a reader checking a
sentence against the figure should find it there. The set is now everything the Results
state about behavior or the partial-activation rate that the sweep records per condition.

What is still missing, and why: the stimulus-, response-signal- and response-window
TRENDS. The sweep records a per-condition trend for the delay window only, so illustrating
those would need the sweep re-run with the contrast computed in every window. The
electroencephalographic measures are absent for a different reason -- the model makes no
quantitative prediction about mu-beta lateralization or the readiness potential, only
about their direction, so there is no range to draw.

Writes results/effect_panel.json.  Usage: python3 effect_panel.py [DUMP.jsonl]
"""
import json, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.targets import build

DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened.jsonl'
OUT = 'results/effect_panel.json'
# cell order (targets.LAB): 0-2 = 0.3 s at 0/12/50% coherence, 3-5 = 2.0 s at 0/12/50%
SHORT, LONG = [0, 1, 2], [3, 4, 5]
C0, C50 = [0, 3], [2, 5]

# Contrast weights, one per effect, applied to the six condition values. Every contrast is
# written so that the PREDICTED sign is positive, so a negative median is a failure of the
# architecture and not a bookkeeping accident, and the names say which difference is taken.
DUR = np.array([1, 1, 1, -1, -1, -1]) / 3.0          # short minus long
COH = np.array([-.5, 0, .5, -.5, 0, .5])             # high minus low coherence
W = {
    # --- behavior
    'acc_coh':   COH,
    'acc_dur':  -DUR,                                 # accuracy is HIGHER after the long
    'acc_int':   np.array([.5, -1, .5, -.5, 1, -.5]),  # duration benefit confined to 12%
    'rt_coh':   -COH,                                 # RT falls as coherence rises
    'rt_dur':    DUR,
    'rt_int':    np.array([1, 0, -1, -1, 0, 1]),      # coherence effect larger after short
    # --- stimulus-onset window
    'stim_coh':  COH,
    'stim_dur':  DUR,
    # --- delay window
    'delay_coh': COH,
    'delay_dur': -DUR,                                # delay rate HIGHER after the long
    'delay_int': np.array([1, 0, -1, -1, 0, 1]),      # coherence effect SMALLER after long
    'trend':     np.ones(6) / 6.0,
    'trend_dur': -DUR,                                # the ramp steepens after the long
    'curv':      np.ones(6) / 6.0,
    # --- response-signal window
    'pmet_coh': -COH,
    'pmet_dur':  DUR,
    'met_dur':   DUR,
}
# Which per-cell quantity each contrast is taken on, and the scale it is reported in.
SRC = {
    'acc_coh': ('acc', 1.0), 'acc_dur': ('acc', 1.0), 'acc_int': ('acc', 1.0),
    'rt_coh': ('rt', 1000.0), 'rt_dur': ('rt', 1000.0), 'rt_int': ('rt', 1000.0),
    'stim_coh': ('srate', 1.0), 'stim_dur': ('srate', 1.0),
    'delay_coh': ('rate', 1.0), 'delay_dur': ('rate', 1.0), 'delay_int': ('rate', 1.0),
    'trend': ('trend', 1.0), 'trend_dur': ('trend', 1.0), 'curv': ('curv', 1.0),
    'pmet_coh': ('pmet', 1000.0), 'pmet_dur': ('pmet', 1000.0), 'met_dur': ('met', 1000.0),
}
# Panel groups, in the order the Results take them. Each entry is (window, [keys]); the
# figure draws one labelled block per window, three panels to a row.
GROUPS = [
    ('Behavior', ['acc_coh', 'acc_dur', 'acc_int', 'rt_coh', 'rt_dur', 'rt_int']),
    ('Stimulus-onset window', ['stim_coh', 'stim_dur']),
    ('Delay window', ['delay_coh', 'delay_dur', 'delay_int',
                      'trend', 'trend_dur', 'curv']),
    ('Response-signal window', ['pmet_coh', 'pmet_dur', 'met_dur']),
]
KEYS = [k for _, ks in GROUPS for k in ks]

NAMES = dict(
    acc_coh='accuracy,\n50% minus 0% coherence',
    acc_dur='accuracy,\nlong minus short',
    acc_int='duration benefit on accuracy,\nat 12% minus at 0% and 50%',
    rt_coh='RT, 0% minus 50% coherence (ms)',
    rt_dur='RT, short minus long (ms)',
    rt_int='coherence effect on RT,\nshort minus long (ms)',
    stim_coh='stimulus-period rate,\n50% minus 0% coherence (/s)',
    stim_dur='stimulus-period rate,\nshort minus long (/s)',
    delay_coh='delay rate,\n50% minus 0% coherence (/s)',
    delay_dur='delay rate,\nlong minus short (/s)',
    delay_int='coherence effect on the delay rate,\nlong minus short (/s)',
    trend='delay-period rate trend (/s$^2$)',
    trend_dur='rate trend, long minus short (/s$^2$)',
    curv='delay-period rate curvature (/s$^3$)',
    pmet_coh='PMET, 0% minus 50% coherence (ms)',
    pmet_dur='PMET, short minus long (ms)',
    met_dur='MET, short minus long (ms)')


def effects(cells):
    """`cells` maps a per-cell quantity name to its six values."""
    return {k: float(np.asarray(W[k], float) @ np.asarray(cells[SRC[k][0]], float)
                     * SRC[k][1]) for k in KEYS}


def effects_se(m):
    """Monte-Carlo standard error of each effect, propagated from the per-cell errors the
    sweep records.

    Without this the sign of an effect is not interpretable. A draw whose true effect is a
    thousandth of a second returns a sign that is a coin toss, and counting it as agreeing
    or disagreeing with the data says nothing about the architecture -- it reports the
    seed. The three-way split (significantly with the observed sign, significantly against
    it, indistinguishable from zero at two standard errors) is what the raw sign fraction
    should have been from the start. It matters most on the delay-period trend, where
    Model 1's value is zero analytically (Appendix, Corollary A.2), so every draw's sign is
    a coin toss and only the three-way split reports that honestly.

    Accuracy gets a binomial error; response time, the two execution times, the trend and
    the curvature use the per-cell standard errors already in the dump; the two rates use
    the Poisson error of the event count behind them, recovered as rate x exposure x
    window.
    """
    g = lambda k: np.asarray(m['raw'][k], float)
    acc, n = g('acc'), g('rt_n')

    # every count-based error is floored at ONE event. A cell in which the simulation
    # produced no premature response at all has an empirical binomial variance of exactly
    # zero, and taking that at face value would declare a difference of a single trial
    # infinitely significant. The floor is the standard small-count device: the
    # uncertainty on an observed zero is the uncertainty on one.
    def bse(p, ntr):
        p = np.maximum(p, 1.0 / ntr)
        return np.sqrt(np.clip(p * (1 - p), 0, None) / ntr)

    def pse(rate, expo):
        return np.sqrt(np.maximum(rate * expo, 1.0)) / np.maximum(expo, 1)

    se = dict(acc=bse(acc, np.maximum(n, 1)),
              rt=g('rt_se'), pmet=g('pmet_se'), met=g('met_se'),
              trend=g('trend_se'), curv=g('curv_se'),
              srate=pse(g('srate'), n * 2.0), rate=pse(g('rate'), n * 2.0))
    out = {}
    for k in KEYS:
        src, scale = SRC[k]
        w = np.asarray(W[k], float)
        out[k] = float(np.sqrt((w ** 2) @ (np.asarray(se[src], float) ** 2))) * scale
    return out


recs = [json.loads(l) for l in open(DUMP)]
E = build()
nsub = E['rt'].shape[0]
obs_cells = [dict(acc=E['acc'][s], rt=E['rt'][s],
                  srate=E['stim'][s].sum(axis=-1) / 2.0,
                  rate=E['dly'][s].sum(axis=-1) / 2.0,
                  trend=E['trend'][s], curv=E['quad'][s],
                  pmet=E['pmet'][s], met=E['met'][s]) for s in range(nsub)]
per_sub = [effects(c) for c in obs_cells]
obs = {k: float(np.nanmean([p[k] for p in per_sub])) for k in KEYS}
obs_se = {k: float(np.nanstd([p[k] for p in per_sub], ddof=1) / np.sqrt(nsub)) for k in KEYS}

out = dict(source=DUMP, names=NAMES, groups=GROUPS, observed=obs, observed_se=obs_se,
           models={})
for tag in ['model1', 'model2']:
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    for need in ('srate', 'trend', 'curv'):
        if adm and need not in adm[0]['raw']:
            raise SystemExit(f"{DUMP} has no `{need}` per cell -- rerun fused_sweep.py")
    rows = [effects({k: np.asarray(m['raw'][src], float) for k, src in
                     (('acc', 'acc'), ('rt', 'rt_mean'), ('srate', 'srate'),
                      ('rate', 'rate'), ('trend', 'trend'), ('curv', 'curv'),
                      ('pmet', 'pmet_mean'), ('met', 'met_mean'))}) for m in adm]
    ses = [effects_se(m) for m in adm]
    out['models'][tag] = dict(n=len(rows))
    for k in KEYS:
        v = np.array([r[k] for r in rows], float)
        e = np.array([x[k] for x in ses], float)
        ok = np.isfinite(v) & np.isfinite(e)
        v, e = v[ok], e[ok]
        # <= and not <, so that a draw whose effect is exactly zero with exactly zero
        # error -- no events in either condition -- counts as uninformative rather than as
        # a disagreement
        zero = np.abs(v) <= 2 * e
        right = (~zero) & (np.sign(v) == np.sign(obs[k]))
        wrong = (~zero) & (np.sign(v) != np.sign(obs[k]))
        out['models'][tag][k] = dict(
            median=float(np.median(v)), p5=float(np.percentile(v, 5)),
            p95=float(np.percentile(v, 95)), min=float(v.min()), max=float(v.max()),
            covers_observed=bool(np.percentile(v, 5) <= obs[k] <= np.percentile(v, 95)),
            sign_agrees=bool(np.sign(np.median(v)) == np.sign(obs[k])),
            frac_right=float(right.mean()), frac_wrong=float(wrong.mean()),
            frac_null=float(zero.mean()), n_tested=int(v.size),
            pct_of_observed=float((v < obs[k]).mean() * 100.0),
            frac_reaching=float(np.mean(v >= obs[k]) if obs[k] > 0
                                else np.mean(v <= obs[k])))

json.dump(out, open(OUT, 'w'), indent=1)
lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
w = max(len(v) for v in lab.values())
print(f"{'effect':{w}} {'observed':>18} {'M1 right/wrong/0':>16} {'M2 right/wrong/0':>16}")
for win, ks in GROUPS:
    print(f"-- {win}")
    for k in ks:
        f = '%+.0f' if SRC[k][1] == 1000.0 else '%+.4f'
        cells = []
        for tag in ['model1', 'model2']:
            r = out['models'][tag][k]
            cells.append(f"{100*r['frac_right']:3.0f}/{100*r['frac_wrong']:3.0f}"
                         f"/{100*r['frac_null']:3.0f}")
        print(f"{lab[k]:{w}} {f % obs[k]:>9} +/-{(f % obs_se[k]).lstrip('+'):>7} "
              f"{cells[0]:>16} {cells[1]:>16}")
print("\ncolumns are the percentage of admissible parameter sets whose effect is, at two "
      "Monte-Carlo\nstandard errors, significantly of the OBSERVED sign / significantly "
      "AGAINST it / not\ndistinguishable from zero")
