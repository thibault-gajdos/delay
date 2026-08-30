"""Ratcliff-style RT-quantile calibration targets (replaces the mean/SD/proportion
targets of the original build_targets.py (removed 2026-08-30) for the
model-vs-behaviour comparison).

Per condition, correct and incorrect responses are treated separately (the classic
quantile-probability decomposition): the 5 Ratcliff quantiles (.1/.3/.5/.7/.9) of RT
partition each response type into 6 bins holding, by construction, 10/20/20/20/20/10 %
of that response type's trials. Premature, sub-0.2-s and no-response-by-4s trials are
dropped from the calibration entirely (not scored, not counted in the denominator):
p(correct) + p(incorrect) = 1 over VALID trials only, exactly as valid was defined in
the original build_targets.py (removed 2026-08-30; see README).

Two-stage aggregation is used where the data support it, one-stage (pooled) where they
do not:
  - p(correct): two-stage -- per-subject accuracy on valid trials, then averaged across
    subjects (identical to the 'acc' target in targets_behav.json).
  - correct-RT quantiles: two-stage -- per-subject quantiles (every subject has 65-151
    correct trials per cell, ample), Vincentized (averaged) across subjects.
  - incorrect-RT quantiles: ONE-STAGE (pooled across all subjects' incorrect trials per
    cell). Per-subject incorrect counts are too sparse for reliable per-subject quantiles
    in the two high-coherence cells (median 3-5 incorrect trials/subject at c=0.50, up to
    16/22 subjects with <10 in a cell) -- per-subject quintiles from that few trials are
    meaningless. Pooling is the standard fallback in this situation.
"""
import numpy as np, json, os

BEH = os.path.join(os.environ.get('DELAY_DATA', '/home/claude/data'), 'data/behavioral')
CELLS = [(0.0, 0.3), (0.12, 0.3), (0.50, 0.3), (0.0, 2.0), (0.12, 2.0), (0.50, 2.0)]
LAB = ['0.3/c00', '0.3/c12', '0.3/c50', '2.0/c00', '2.0/c12', '2.0/c50']
NS = 22
Q = [0.1, 0.3, 0.5, 0.7, 0.9]
MIN_SUBJ_Q = 20   # minimum per-subject trials to include that subject in Vincentizing

cells_out = []
for k, (c, D) in enumerate(CELLS):
    subj_acc = []
    subj_cor_rt = []          # list of per-subject correct-RT arrays (for Vincentizing)
    pooled_incor_rt = []
    n_valid_tot = 0; n_cor_tot = 0; n_incor_tot = 0
    for s in range(NS):
        b = np.load(f'{BEH}/s{s+1}.npy', allow_pickle=True)
        coh = b[:, 5].astype(float); dur = b[:, 6].astype(float).astype(int)
        rm = b[:, 10].astype(float).astype(int); ac = b[:, 11].astype(float).astype(int)
        rt = b[:, 12].astype(float); tg = b[:, 14].astype(float).astype(int)
        sel = (coh == c) & (dur == (18 if D == 0.3 else 120))
        responded = sel & (rm == 1) & (ac != -1)
        valid = responded & (rt >= 0.2) & (rt <= 4.0)
        cor = valid & (ac == 1); incor = valid & (ac == 0)
        n_valid_tot += valid.sum(); n_cor_tot += cor.sum(); n_incor_tot += incor.sum()
        subj_acc.append(ac[valid].mean())
        if cor.sum() >= MIN_SUBJ_Q:
            subj_cor_rt.append(np.percentile(rt[cor], [q * 100 for q in Q]))
        pooled_incor_rt.append(rt[incor])

    p_correct = float(np.mean(subj_acc))
    cor_q = np.mean(np.vstack(subj_cor_rt), axis=0)          # Vincentized
    incor_all = np.concatenate(pooled_incor_rt)
    incor_q = np.percentile(incor_all, [q * 100 for q in Q]) if incor_all.size >= 20 \
        else np.full(5, np.nan)

    cells_out.append(dict(
        coh=c, dur=D, label=LAB[k],
        n_valid=int(n_valid_tot), n_correct=int(n_cor_tot), n_incorrect=int(n_incor_tot),
        n_subj_vincent_correct=len(subj_cor_rt),
        p_correct=p_correct, p_incorrect=1.0 - p_correct,
        correct_quantiles=cor_q.tolist(), incorrect_quantiles=incor_q.tolist(),
    ))

out = dict(cells=cells_out, quantile_levels=Q, nsub=NS,
           note="premature / <0.2s / no-response-by-4s trials excluded entirely; "
                "p_correct+p_incorrect=1 over valid trials only.")
json.dump(out, open('results/targets_q.json', 'w'), indent=1)

print(f"{'cell':>8} {'p_cor':>6} {'n_val':>6} {'n_cor':>6} {'n_inc':>6} {'#subjVinc':>9}  "
      f"correct Q(.1/.3/.5/.7/.9)  |  incorrect Q")
for cc in cells_out:
    cq = " ".join(f"{x:.3f}" for x in cc['correct_quantiles'])
    iq = " ".join(f"{x:.3f}" for x in cc['incorrect_quantiles'])
    print(f"{cc['label']:>8} {cc['p_correct']:6.3f} {cc['n_valid']:6d} {cc['n_correct']:6d} "
          f"{cc['n_incorrect']:6d} {cc['n_subj_vincent_correct']:9d}  {cq}  |  {iq}")
