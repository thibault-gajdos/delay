"""
effect_panel.py -- the effect sizes each architecture can produce, on the predictions that
actually discriminate, against the observed value.

This is what the parameter-space evaluation is for. The predictions of Appendices A-E are
ordinal and hold at every admissible setting; what a derivation cannot supply is the scale.
For each effect below, this reports the distribution of that effect over a model's
admissible draws -- median and 5th-95th percentile -- next to the observed value and its
between-participant standard error.

The effects are chosen to be the ones that carry an argument, which means the contrasts and
the interactions rather than the main effects. A main effect that both architectures produce
with the same sign at every parameter setting tells the reader nothing he did not already
know from Appendix A; what is worth reporting is (i) the interactions, which is where the
architectures' shared machinery makes a specific and falsifiable claim, (ii) the one
prediction outside the delay period on which the two models disagree, and (iii) the
dissociation between PMET and MET, which is the claim that gives the response-time
decomposition its meaning.

  A.2   duration benefit on accuracy, at 12% coherence minus at 0% and 50%
        The duration benefit is predicted to vanish at both ends of the coherence range
        (accuracy is exactly 1/2 at c = 0 whatever the duration, and tends to 1 at high c
        for both durations), so it must be confined to intermediate coherence. A positive
        value is the prediction.
  RT    duration effect, and coherence effect (ms)
        The two headline behavioural effects, kept because their SIZE is the thing a
        derivation cannot give.
  B.2   duration effect on the stimulus-period PMA rate (short minus long, /s)
        THE discriminating prediction outside the delay: Model 2 predicts a higher and more
        steeply rising stimulus-locked rate after the SHORT stimulus, Model 1 the opposite.
        The sign is what matters here, not the magnitude.
  D     duration effect on PMET, and on MET (ms)
        Reported together because the claim is the dissociation: the conditions move the
        premotor component and leave the motor component alone. An architecture that moves
        both, or neither, has not reproduced the result whatever either effect's size.
  C.3   coherence x duration interaction on the delay-period PMA rate (/s)
        Saturation of the maintained level predicts the coherence effect to SHRINK after a
        long stimulus. The data show it growing. This is the paper's acknowledged
        disagreement, and it belongs in the panel for that reason.
  C.2   delay-period log-rate slope (/s)
        The quantity on which the two architectures separate.

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


def effects(acc, rt, srate, rate, pmet, met, prem, slope):
    """Every directional prediction the Results test, in the sign the proposition predicts.

    Each entry is defined so that the PREDICTED sign is positive, so a negative median is a
    failure of the architecture and not a bookkeeping accident. The names say which
    difference is taken, because a duration effect can be written either way round and the
    manuscript's propositions are statements about a direction.
    """
    return {
        # A.2: the duration benefit on accuracy is confined to intermediate coherence
        'acc_int': float((acc[4] - acc[1]) - 0.5 * ((acc[3] - acc[0]) + (acc[5] - acc[2]))),
        # A.3: response time is shorter after the long stimulus and at high coherence
        'rt_dur': float(np.nanmean(rt[SHORT]) - np.nanmean(rt[LONG])) * 1000.0,
        # A.4: the coherence effect on response time is larger after the short stimulus
        'rt_int': float((rt[0] - rt[2]) - (rt[3] - rt[5])) * 1000.0,
        # B.1: the stimulus-period rate is higher at high coherence (both models)
        'stim_coh': float(np.nanmean(srate[C50]) - np.nanmean(srate[C0])),
        # B.2: the stimulus-period rate is higher after the SHORT stimulus under Model 2 and
        # after the long one under Model 1 -- the one prediction outside the delay on which
        # they are supposed to disagree
        'stim_dur': float(np.nanmean(srate[SHORT]) - np.nanmean(srate[LONG])),
        # D.2: the conditions move the premotor component and leave the motor one alone
        'pmet_dur': float(np.nanmean(pmet[SHORT]) - np.nanmean(pmet[LONG])) * 1000.0,
        'met_dur': float(np.nanmean(met[SHORT]) - np.nanmean(met[LONG])) * 1000.0,
        # C.3: saturation makes the coherence effect on the delay rate SMALLER after the long
        # stimulus, so the predicted sign of this interaction is negative; it is written
        # long-minus-short so that, as everywhere else here, the observed sign is positive
        'delay_int': float((rate[5] - rate[3]) - (rate[2] - rate[0])),
        # C.2: the delay-period log-rate slope
        'slope': float(np.nanmean(slope)),
    }


# Panel order fixes the panel letters A-I that the figure draws and the Results cite.
# The old "A.2 / B.1 / C.3 ..." prefixes named statements in the five-appendix scheme that
# was collapsed into a single appendix on 2026-08-27; they were stale and are gone.
NAMES = dict(
    acc_int='duration benefit on accuracy,\nat 12% minus at 0% and 50%',
    rt_dur='RT, short minus long (ms)',
    rt_int='coherence effect on RT,\nshort minus long (ms)',
    stim_coh='stimulus-period rate,\n50% minus 0% coherence (/s)',
    stim_dur='stimulus-period rate,\nshort minus long (/s)',
    pmet_dur='PMET, short minus long (ms)',
    met_dur='MET, short minus long (ms)',
    delay_int='coherence effect on the delay rate,\nlong minus short (/s)',
    slope='delay-period log-rate slope (/s)')

NTR = 4000.0


def effects_se(m):
    """Monte-Carlo standard error of each effect, propagated from the per-cell errors the
    sweep records.

    Without this the sign of an effect is not interpretable. A draw whose true effect is a
    thousandth of a second returns a sign that is a coin toss, and counting it as agreeing or
    disagreeing with the data says nothing about the architecture -- it reports the seed. The
    three-way split below (significantly with the observed sign, significantly against it,
    indistinguishable from zero at two standard errors) is what the raw sign fraction should
    have been from the start: on the delay-period slope it turns Model 1's "31% with the
    observed sign" into "3% significantly positive, 39% significantly negative, 58% not
    distinguishable from zero", which is a different and more honest statement.

    Accuracy gets a binomial error; response time, the two execution times and the
    log-rate slope use the per-cell standard errors already in the dump; the two rates use
    the Poisson error of the event count behind them, recovered as rate x exposure x
    window.
    """
    g = lambda k: np.asarray(m['raw'][k], float)
    acc, rt, n = g('acc'), g('rt_mean'), g('rt_n')
    # every count-based error is floored at ONE event. A cell in which the simulation
    # produced no premature response at all has an empirical binomial variance of exactly
    # zero, and taking that at face value would declare a difference of a single trial
    # infinitely significant -- which is how a draw with no premature responses anywhere
    # came to be counted as contradicting the data. The floor is the standard
    # small-count device: the uncertainty on an observed zero is the uncertainty on one.
    def bse(p, ntr):
        p = np.maximum(p, 1.0 / ntr)
        return np.sqrt(np.clip(p * (1 - p), 0, None) / ntr)

    def pse(rate, expo):
        # Poisson error of the event count behind a rate, recovered as rate x exposure
        return np.sqrt(np.maximum(rate * expo, 1.0)) / np.maximum(expo, 1)

    accse = bse(acc, np.maximum(n, 1))
    srse = pse(g('srate'), n * 2.0)
    rse2 = pse(g('rate'), n * 2.0)

    def lin(se, w, scale=1.0):
        w = np.asarray(w, float)
        return float(np.sqrt((w ** 2) @ (se ** 2))) * scale

    return dict(
        acc_int=lin(accse, [.5, -1, .5, -.5, 1, -.5]),
        rt_dur=lin(g('rt_se'), np.array([1, 1, 1, -1, -1, -1]) / 3.0, 1000.0),
        rt_int=lin(g('rt_se'), [1, 0, -1, -1, 0, 1], 1000.0),
        stim_coh=lin(srse, [-.5, 0, .5, -.5, 0, .5]),
        stim_dur=lin(srse, np.array([1, 1, 1, -1, -1, -1]) / 3.0),
        pmet_dur=lin(g('pmet_se'), np.array([1, 1, 1, -1, -1, -1]) / 3.0, 1000.0),
        met_dur=lin(g('met_se'), np.array([1, 1, 1, -1, -1, -1]) / 3.0, 1000.0),
        delay_int=lin(rse2, [1, 0, -1, -1, 0, 1]),
        slope=lin(g('slope_se'), np.ones(6) / 6.0))


recs = [json.loads(l) for l in open(DUMP)]
E = build()
nsub = E['rt'].shape[0]
stim_cells = E['stim'].sum(axis=-1) / 2.0        # (nsub, 6) rate over the 0-2 s window
rate_cells = E['dly'].sum(axis=-1) / 2.0

per_sub = [effects(E['acc'][s], E['rt'][s], stim_cells[s], rate_cells[s], E['pmet'][s],
                   E['met'][s], E['prem'][s], E['slope'][s]) for s in range(nsub)]
obs = {k: float(np.nanmean([p[k] for p in per_sub])) for k in NAMES}
obs_se = {k: float(np.nanstd([p[k] for p in per_sub], ddof=1) / np.sqrt(nsub)) for k in NAMES}

out = dict(source=DUMP, names=NAMES, observed=obs, observed_se=obs_se, models={})
for tag in ['model1', 'model2']:
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    if adm and 'srate' not in adm[0]['raw']:
        raise SystemExit(f"{DUMP} predates the srate/pmet fields -- rerun fused_sweep.py")
    rows = [effects(*[np.asarray(m['raw'][k], float) for k in
                      ('acc', 'rt_mean', 'srate', 'rate', 'pmet_mean', 'met_mean', 'prem',
                       'slope')]) for m in adm]
    ses = [effects_se(m) for m in adm]
    out['models'][tag] = dict(n=len(rows))
    for k in NAMES:
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
            # the fraction of admissible draws whose effect has the observed SIGN, and the
            # fraction reaching the observed magnitude. For B.2 the sign is what
            # discriminates -- Proposition B.2 says the two models predict opposite signs --
            # so a range that fails to cover the observed value can still be the decisive
            # comparison, and reporting only the interval would hide that.
            frac_right=float(right.mean()), frac_wrong=float(wrong.mean()),
            frac_null=float(zero.mean()), n_tested=int(v.size),
            frac_reaching=float(np.mean(v >= obs[k]) if obs[k] > 0 else np.mean(v <= obs[k])))

json.dump(out, open(OUT, 'w'), indent=1)
lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
w = max(len(v) for v in lab.values())
print(f"{'effect':{w}} {'observed':>16} {'M1 right/wrong/0':>16} {'M2 right/wrong/0':>16}")
for k in NAMES:
    f = '%+.0f' if k in ('rt_dur', 'rt_int', 'pmet_dur', 'met_dur') else '%+.4f'
    cells = []
    for tag in ['model1', 'model2']:
        r = out['models'][tag][k]
        cells.append(f"{100*r['frac_right']:3.0f}/{100*r['frac_wrong']:3.0f}/{100*r['frac_null']:3.0f}")
    print(f"{lab[k]:{w}} {f % obs[k]} +/-{(f % obs_se[k]).lstrip('+')} "
          f"{cells[0]:>16} {cells[1]:>16}")
print("\ncolumns are the percentage of admissible parameter sets whose effect is, at two "
      "Monte-Carlo\nstandard errors, significantly of the OBSERVED sign / significantly "
      "AGAINST it / not\ndistinguishable from zero")
