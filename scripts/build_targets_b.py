"""build_targets_b.py -- targets for the ALL-DATA fit (fit B).

Fit A (build_targets_q.py -> targets_q.json, calib_behav_q.py) is unchanged and remains
the manuscript's primary fit: behaviour only, with the delay-period log-rate slope held
out as an out-of-sample check.  Fit B is additional and sits at a different level: it asks
whether one parameter vector can reproduce behaviour AND the electromyographic measures
together.  It changes no prediction -- the predictions are derived analytically and tested
against the data, not against either fit.

Families per cell (35 categories, against fit A's 13).  Every family is a multinomial
scored by the same deviance, so they add on one scale with no free weighting constant:

  premature                  1 category   (as in fit A)
  PMET quantiles x {cor,inc} 6 + 6        (replaces fit A's RT quantiles -- RT is the
                                           exact sum PMET + MET, so keeping all three
                                           would count the same information twice)
  MET quantiles  x {cor,inc} 6 + 6
  delay PMA count            4            {0, 1, 2, >=3} partial activations per trial
  delay PMA shape            5            onset time over 5 bins of 0.4 s, given >=1
  stimulus / RS / response   2 + 2 + 2    P(>=1 partial activation) in each window

Why counts rather than a Poisson deviance.  Partial-activation counts are massively
overdispersed between participants -- variance/mean ~= 100 in the delay window, against 1
for a Poisson -- so a Poisson deviance would give the EMG block roughly a hundred times
its real information.  Categorising the count absorbs the overdispersion with no
distributional assumption and reuses the multinomial machinery unchanged.

Why {0,1,2,>=3} only for the delay.  Measured per-trial count distributions:
    stimulus  mean 0.47   P(0)=83.1%  P(1)= 8.4%  P(2)=3.0%  P(>=3)= 5.4%
    delay     mean 1.42   P(0)=65.6%  P(1)=11.1%  P(2)=5.4%  P(>=3)=17.9%
    RS        mean 0.27   P(0)=80.0%  P(1)=14.3%  P(2)=5.0%  P(>=3)= 0.8%
    response  mean 0.02   P(0)=97.5%  P(1)= 2.4%  P(2)=0.0%  P(>=3)= 0.0%
At a pseudo-sample size of 22 a 5% category carries 1.1 pseudo-observations, so only the
delay window supports four categories; the others collapse to P(>=1).  The response
window is nearly empty (97.5% zeros) and will contribute almost nothing -- it is included
for completeness, not because it is expected to constrain anything.

Why 5 bins for the delay shape.  At pseudo-N = 22 and 1.42 activations per trial, the
original 20 bins of 0.1 s give 1.56 expected per bin, too sparse for a stable G^2; 5 bins
of 0.4 s give 6.2 and still resolve the ramp cleanly (observed shape, pooled:
0.094 0.126 0.182 0.254 0.344 against 0.200 for a flat profile).

Aggregation follows fit A: two stages (per participant, then averaged) wherever the data
support it, one stage (pooled) where they do not.  Correct-response quantiles and all the
PMA families are two-stage; incorrect-response quantiles are pooled, because at the higher
coherences several participants have too few errors for a per-participant quantile.

NOTE ON DENOMINATORS.  empirical.trials() already restricts to valid trials and
additionally drops trials with missing triggers or no identifiable response burst, so the
EMG families rest on a slightly smaller trial set than the behavioural ones.  Each family
is internally consistent; they are not forced onto a common denominator.

Writes results/targets_b.json.
"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.empirical import trials

Q = [0.1, 0.3, 0.5, 0.7, 0.9]
CELLS = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]   # (coh idx, dur idx)
LAB = {(0, 0): '0.3/c00', (1, 0): '0.3/c12', (2, 0): '0.3/c50',
       (0, 1): '2.0/c00', (1, 1): '2.0/c12', (2, 1): '2.0/c50'}
COHV = {0: 0.0, 1: 0.12, 2: 0.50}
DURV = {0: 0.3, 1: 2.0}
NSUB = 22
MIN_Q = 20            # per-participant trials needed to enter the Vincentized average
N_SHAPE = 5           # delay-window shape bins
DELAY_W = (1.0, 3.0)  # relative to stimulus offset


def windows(tr):
    """per-trial partial-activation counts in the four analysis windows, and the delay
    onset times measured from the start of the delay window"""
    pt = np.asarray(tr['pma_t'], float)
    off, rs, rp = tr['t_off'], tr['t_rs'], tr['t_resp']
    n = lambda a, b: int(((pt >= a) & (pt < b)).sum())
    dly = pt[(pt >= off + DELAY_W[0]) & (pt < off + DELAY_W[1])] - (off + DELAY_W[0])
    return dict(stim=n(0.0, 2.0), delay=dly.size, rs=n(rs, rs + 0.2),
                resp=n(rp - 0.2, rp), delay_onsets=dly)


# ---- gather -----------------------------------------------------------------
per = {c: [] for c in CELLS}          # per-subject records
for sid in range(1, NSUB + 1):
    bysub = {c: dict(pmet_c=[], pmet_i=[], met_c=[], met_i=[],
                     cnt=[], stim=[], rs=[], resp=[], onsets=[]) for c in CELLS}
    for tr in trials(sid):
        key = (tr['coh'], tr['dur'])
        if key not in bysub:
            continue
        w = windows(tr)
        r = bysub[key]
        pmet = tr['t_burst'] - tr['t_rs']
        met = tr['t_resp'] - tr['t_burst']
        if np.isfinite(pmet) and np.isfinite(met) and pmet > 0 and met > 0:
            (r['pmet_c'] if tr['acc'] == 1 else r['pmet_i']).append(pmet)
            (r['met_c'] if tr['acc'] == 1 else r['met_i']).append(met)
        r['cnt'].append(w['delay'])
        r['stim'].append(w['stim'] >= 1)
        r['rs'].append(w['rs'] >= 1)
        r['resp'].append(w['resp'] >= 1)
        r['onsets'].extend(w['delay_onsets'])
    for c in CELLS:
        per[c].append(bysub[c])

# ---- reduce -----------------------------------------------------------------
cells_out = []
for c in CELLS:
    R = per[c]
    def vincentize(field):
        rows = [np.percentile(r[field], [q * 100 for q in Q])
                for r in R if len(r[field]) >= MIN_Q]
        return np.mean(np.vstack(rows), axis=0), len(rows)
    def pooled(field):
        allv = np.concatenate([np.asarray(r[field], float) for r in R]) \
            if any(len(r[field]) for r in R) else np.array([])
        return (np.percentile(allv, [q * 100 for q in Q]) if allv.size >= MIN_Q
                else np.full(len(Q), np.nan)), allv.size

    pmet_c, n_pc = vincentize('pmet_c')
    met_c, n_mc = vincentize('met_c')
    pmet_i, n_pi = pooled('pmet_i')
    met_i, n_mi = pooled('met_i')

    # accuracy split, two-stage
    p_cor = float(np.mean([len(r['pmet_c']) / max(len(r['pmet_c']) + len(r['pmet_i']), 1)
                           for r in R]))

    # delay count distribution {0,1,2,>=3}, two-stage
    def cnt_dist(r):
        a = np.asarray(r['cnt'], int)
        return np.array([np.mean(a == 0), np.mean(a == 1), np.mean(a == 2), np.mean(a >= 3)])
    count_p = np.mean(np.vstack([cnt_dist(r) for r in R if len(r['cnt'])]), axis=0)

    # delay shape over N_SHAPE bins, two-stage over participants with any activation
    rows = []
    for r in R:
        o = np.asarray(r['onsets'], float)
        if o.size >= 10:
            h, _ = np.histogram(o, bins=N_SHAPE, range=(0.0, DELAY_W[1] - DELAY_W[0]))
            rows.append(h / h.sum())
    shape_p = np.mean(np.vstack(rows), axis=0)

    p_any = {k: float(np.mean([np.mean(r[k]) for r in R if len(r[k])]))
             for k in ('stim', 'rs', 'resp')}

    cells_out.append(dict(
        label=LAB[c], coh=COHV[c[0]], dur=DURV[c[1]],
        p_correct=p_cor, p_incorrect=1.0 - p_cor,
        pmet_quantiles_correct=pmet_c.tolist(), pmet_quantiles_incorrect=pmet_i.tolist(),
        met_quantiles_correct=met_c.tolist(), met_quantiles_incorrect=met_i.tolist(),
        delay_count_p=count_p.tolist(), delay_shape_p=shape_p.tolist(),
        p_any_stim=p_any['stim'], p_any_rs=p_any['rs'], p_any_resp=p_any['resp'],
        n_subj_pmet_correct=n_pc, n_pooled_pmet_incorrect=int(n_pi),
        n_trials=int(sum(len(r['cnt']) for r in R)),
    ))

# Incorrect-response chronometry is pooled across ALL cells, not per cell.  At the two
# 50%-coherence cells errors are ~5% of trials in the data and rarer still in the model,
# so a per-cell 6-bin conditional distribution cannot be estimated on either side; empty
# simulated bins would then hit the G^2 floor and charge ~30 per bin, penalising the
# error RATE a second time when the response-type family already scores it.  Pooling is
# the same sparsity logic already applied across participants, extended one level.
allsub_pmet_i, allsub_met_i = [], []
for c in CELLS:
    for r in per[c]:
        allsub_pmet_i.extend(r['pmet_i']); allsub_met_i.extend(r['met_i'])
pool_pmet_i = np.percentile(np.asarray(allsub_pmet_i, float), [q*100 for q in Q])
pool_met_i = np.percentile(np.asarray(allsub_met_i, float), [q*100 for q in Q])

out = dict(cells=cells_out, quantile_levels=Q, nsub=NSUB, n_shape_bins=N_SHAPE,
           delay_window=list(DELAY_W),
           pmet_quantiles_incorrect_pooled=pool_pmet_i.tolist(),
           met_quantiles_incorrect_pooled=pool_met_i.tolist(),
           n_incorrect_pooled=len(allsub_pmet_i),
           note='fit B targets; see build_targets_b.py docstring')
json.dump(out, open('results/targets_b.json', 'w'), indent=1)

print(f"{'cell':>9} {'ntr':>5} {'p_cor':>6} {'PMET .5':>8} {'MET .5':>7}  "
      f"{'count {0,1,2,3+}':>26}  {'shape (5 bins)':>32}")
for cc in cells_out:
    cp = " ".join(f"{x:.3f}" for x in cc['delay_count_p'])
    sp = " ".join(f"{x:.3f}" for x in cc['delay_shape_p'])
    print(f"{cc['label']:>9} {cc['n_trials']:5d} {cc['p_correct']:6.3f} "
          f"{cc['pmet_quantiles_correct'][2]:8.3f} {cc['met_quantiles_correct'][2]:7.3f}  "
          f"{cp:>26}  {sp:>32}")
print(f"\n{'cell':>9} {'P(>=1) stim':>12} {'RS':>8} {'resp':>8}")
for cc in cells_out:
    print(f"{cc['label']:>9} {cc['p_any_stim']:12.3f} {cc['p_any_rs']:8.3f} "
          f"{cc['p_any_resp']:8.3f}")
