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


def effects(acc, rt, srate, pmet, met, rate, slope):
    """The seven effects, from per-cell arrays. Durations in seconds, rates per second."""
    ben = [acc[3 + i] - acc[i] for i in range(3)]         # duration benefit per coherence
    return dict(
        acc_int=float(ben[1] - 0.5 * (ben[0] + ben[2])),
        rt_dur=float(np.nanmean(rt[SHORT]) - np.nanmean(rt[LONG])) * 1000.0,
        rt_coh=float(np.nanmean(rt[[0, 3]]) - np.nanmean(rt[[2, 5]])) * 1000.0,
        stim_dur=float(np.nanmean(srate[SHORT]) - np.nanmean(srate[LONG])),
        pmet_dur=float(np.nanmean(pmet[SHORT]) - np.nanmean(pmet[LONG])) * 1000.0,
        met_dur=float(np.nanmean(met[SHORT]) - np.nanmean(met[LONG])) * 1000.0,
        delay_int=float((rate[5] - rate[3]) - (rate[2] - rate[0])),
        slope=float(np.nanmean(slope)))


NAMES = dict(
    acc_int='duration benefit on accuracy,\nconfined to 12% coherence',
    rt_dur='duration effect on RT (ms)',
    rt_coh='coherence effect on RT (ms)',
    stim_dur='duration effect on the\nstimulus-period PMA rate (/s)',
    pmet_dur='duration effect on PMET (ms)',
    met_dur='duration effect on MET (ms)',
    delay_int='coherence x duration on the\ndelay-period PMA rate (/s)',
    slope='delay-period log-rate slope (/s)')

recs = [json.loads(l) for l in open(DUMP)]
E = build()
nsub = E['rt'].shape[0]
stim_cells = E['stim'].sum(axis=-1) / 2.0        # (nsub, 6) rate over the 0-2 s window
rate_cells = E['dly'].sum(axis=-1) / 2.0

per_sub = [effects(E['acc'][s], E['rt'][s], stim_cells[s], E['pmet'][s], E['met'][s],
                   rate_cells[s], E['slope'][s]) for s in range(nsub)]
obs = {k: float(np.nanmean([p[k] for p in per_sub])) for k in NAMES}
obs_se = {k: float(np.nanstd([p[k] for p in per_sub], ddof=1) / np.sqrt(nsub)) for k in NAMES}

out = dict(source=DUMP, names=NAMES, observed=obs, observed_se=obs_se, models={})
for tag in ['model1', 'model2']:
    adm = [r[tag] for r in recs if r[tag]['admissible']]
    if adm and 'srate' not in adm[0]['raw']:
        raise SystemExit(f"{DUMP} predates the srate/pmet fields -- rerun fused_sweep.py")
    rows = [effects(*[np.asarray(m['raw'][k], float) for k in
                      ('acc', 'rt_mean', 'srate', 'pmet_mean', 'met_mean', 'rate', 'slope')])
            for m in adm]
    out['models'][tag] = dict(n=len(rows))
    for k in NAMES:
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
        out['models'][tag][k] = dict(
            median=float(np.median(v)), p5=float(np.percentile(v, 5)),
            p95=float(np.percentile(v, 95)), min=float(v.min()), max=float(v.max()),
            covers_observed=bool(np.percentile(v, 5) <= obs[k] <= np.percentile(v, 95)),
            sign_agrees=bool(np.sign(np.median(v)) == np.sign(obs[k])))

json.dump(out, open(OUT, 'w'), indent=1)
lab = {k: v.replace('\n', ' ') for k, v in NAMES.items()}
w = max(len(v) for v in lab.values())
print(f"{'effect':{w}} {'observed':>16} {'Model 1  5-95%':>24} {'Model 2  5-95%':>24}")
for k in NAMES:
    f = '%+.3f' if k in ('acc_int', 'stim_dur', 'delay_int', 'slope') else '%+.0f'
    cells = []
    for tag in ['model1', 'model2']:
        r = out['models'][tag][k]
        cells.append(f"{f % r['p5']} to {f % r['p95']}{'' if r['covers_observed'] else '*'}")
    print(f"{lab[k]:{w}} {f % obs[k]} +/-{(f % obs_se[k]).lstrip('+')} "
          f"{cells[0]:>24} {cells[1]:>24}")
print("\n* = the model's 5th-95th percentile does not cover the observed value")
