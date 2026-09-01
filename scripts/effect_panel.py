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
        # A.5: premature responses are more frequent after the long stimulus
        'prem_dur': float(np.nanmean(prem[LONG]) - np.nanmean(prem[SHORT])),
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


NAMES = dict(
    acc_int='A.2  duration benefit on accuracy,\nat 12% minus at 0% and 50%',
    rt_dur='A.3  RT, short minus long (ms)',
    rt_int='A.4  coherence effect on RT,\nshort minus long (ms)',
    prem_dur='A.5  premature responses,\nlong minus short',
    stim_coh='B.1  stimulus-period rate,\n50% minus 0% coherence (/s)',
    stim_dur='B.2  stimulus-period rate,\nshort minus long (/s)',
    pmet_dur='D.2  PMET, short minus long (ms)',
    met_dur='D.2  MET, short minus long (ms)',
    delay_int='C.3  coherence effect on the delay rate,\nlong minus short (/s)',
    slope='C.2  delay-period log-rate slope (/s)')

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
    out['models'][tag] = dict(n=len(rows))
    for k in NAMES:
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
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
            # computed over the sets whose effect is non-zero: a draw with exactly zero
            # premature responses in both conditions carries no directional information,
            # and counting it as disagreeing would understate the agreement
            frac_correct_sign=(float(np.mean(np.sign(v[v != 0]) == np.sign(obs[k])))
                               if np.any(v != 0) else float('nan')),
            n_nonzero=int(np.sum(v != 0)),
            frac_reaching=float(np.mean(v >= obs[k]) if obs[k] > 0 else np.mean(v <= obs[k])))

json.dump(out, open(OUT, 'w'), indent=1)
lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
w = max(len(v) for v in lab.values())
print(f"{'effect':{w}} {'observed':>16} {'Model 1  5-95%':>24} {'Model 2  5-95%':>24}")
for k in NAMES:
    f = '%+.0f' if k in ('rt_dur', 'rt_int', 'pmet_dur', 'met_dur') else '%+.4f'
    cells = []
    for tag in ['model1', 'model2']:
        r = out['models'][tag][k]
        cells.append(f"{f % r['p5']} to {f % r['p95']}{'' if r['covers_observed'] else '*'}")
    print(f"{lab[k]:{w}} {f % obs[k]} +/-{(f % obs_se[k]).lstrip('+')} "
          f"{cells[0]:>24} {cells[1]:>24}")
print("\n* = the model's 5th-95th percentile does not cover the observed value")
