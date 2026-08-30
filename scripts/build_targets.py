"""Behavioural calibration targets, two-stage aggregation.

Stage (i): within participant, within condition -> one number per participant x cell.
Stage (ii): across participants -> the target, plus its between-participant SE.

Trial classification, from the raw behavioural files (col 10 = response moment,
col 11 = accuracy, col 12 = RT from RS onset, col 14 = response trigger):
    pre-RS      col10 == 0                      (response before the response signal)
    no-response col14 == 250                    (no response within the 4-s deadline)
    fast        valid response with RT < 0.2 s
    valid       everything else; accuracy and RT use these trials ONLY
"""
import numpy as np, json, os
BEH = os.path.join(os.environ.get('DELAY_DATA', os.path.expanduser('~/mnt/delay/data/86cbw-osfstorage-archive')), 'data/behavioral')
CELLS = [(0.0, 0.3), (0.12, 0.3), (0.50, 0.3), (0.0, 2.0), (0.12, 2.0), (0.50, 2.0)]
LAB = ['0.3/c00', '0.3/c12', '0.3/c50', '2.0/c00', '2.0/c12', '2.0/c50']
NS = 22

blocks = {b: np.zeros((NS, 6)) for b in ('acc', 'rt', 'rt_iqr', 'prem', 'fast', 'noresp')}
ntr = np.zeros((NS, 6))
for s in range(NS):
    b = np.load(f'{BEH}/s{s+1}.npy', allow_pickle=True)
    coh = b[:, 5].astype(float); dur = b[:, 6].astype(float).astype(int)
    rm = b[:, 10].astype(float).astype(int); ac = b[:, 11].astype(float).astype(int)
    rt = b[:, 12].astype(float); tg = b[:, 14].astype(float).astype(int)
    for k, (c, D) in enumerate(CELLS):
        sel = (coh == c) & (dur == (18 if D == 0.3 else 120))
        ntr[s, k] = sel.sum()
        blocks['prem'][s, k]   = (rm[sel] == 0).mean()
        blocks['noresp'][s, k] = (tg[sel] == 250).mean()
        responded = sel & (rm == 1) & (ac != -1)
        blocks['fast'][s, k]   = (responded & (rt < 0.2)).sum() / sel.sum()
        v = responded & (rt >= 0.2) & (rt <= 4.0)
        blocks['acc'][s, k]    = ac[v].mean()
        blocks['rt'][s, k]     = rt[v].mean()
        blocks['rt_iqr'][s, k] = np.subtract(*np.percentile(rt[v], [75, 25]))

out = {'cells': CELLS, 'lab': LAB, 'nsub': NS, 'ntr': ntr.tolist()}
for name, A in blocks.items():
    out[name] = {'per_subject': A.tolist(),
                 'target': A.mean(0).tolist(),
                 'se': (A.std(0, ddof=1) / np.sqrt(NS)).tolist()}
json.dump(out, open('results/targets_behav.json', 'w'))

print(f"{'':>10}" + "".join(f"{l:>10}" for l in LAB))
for name in ('acc', 'rt', 'rt_iqr', 'prem', 'fast', 'noresp'):
    t = np.array(out[name]['target']); e = np.array(out[name]['se'])
    print(f"{name:>10}" + "".join(f"{x:10.4f}" for x in t))
    print(f"{'  (SE)':>10}" + "".join(f"{x:10.4f}" for x in e))
print("\ntrials per cell (mean):", ntr.mean(0).round(1))
