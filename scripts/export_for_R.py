"""
export_for_R.py -- the bridge between the Python side and the R side.

Python owns the model, the simulation and the extraction of events from the
electromyographic files; R owns every statistical test.  This script writes the
one file R expects:

    pma_counts.csv      empirical counts, one row per participant, cell, period
                        and bin, with the raw count and the number of trials
                        contributing to that bin.  Needed for the Poisson
                        analyses, which require counts and exposures rather than
                        rates.

Until 2026-09-05 it also wrote sim_counts_m{1,2}.csv, the same quantities from each
model simulated at the representative parameter vector.  Nothing read them: R analyses
the empirical counts only, and what the models predict is reported as effect-size ranges
over the parameter space (effect_sizes.py, effect_panel.py).  Both the files and the
function that built them are gone.

Exposure is not always the number of valid trials.  A partial activation is by
definition an activation preceding the response burst, so once that burst has
begun a trial can no longer contribute, and bins after it are undefined rather
than empty.  This matters only in the response-signal and response windows,
where a few per cent of trials have already produced their burst; it is handled
by counting, for each bin, the trials whose burst had not yet begun.
"""
import numpy as np, sys, csv, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.empirical import trials

NSUB = 22                      # participants in the experiment
PERIODS = ['stimulus', 'delay', 'RS', 'response']


def empirical_rows():
    rows = []
    for s in range(NSUB):
        cnt = {p: np.zeros((6, 20 if p in ('stimulus', 'delay') else 2)) for p in PERIODS}
        expo = {p: np.zeros((6, 20 if p in ('stimulus', 'delay') else 2)) for p in PERIODS}
        for tr in trials(s + 1):
            k = tr['dur'] * 3 + tr['coh']
            v = tr['pma_t']
            u = v[(v >= 0) & (v < 2.0) & (v < tr['t_burst'])]
            cnt['stimulus'][k] += np.histogram(u, bins=20, range=(0, 2))[0]
            expo['stimulus'][k] += 1
            w = v[(v >= tr['t_off'] + 1) & (v < tr['t_off'] + 3)] - tr['t_off'] - 1
            cnt['delay'][k] += np.histogram(w, bins=20, range=(0, 2))[0]
            expo['delay'][k] += 1
            q = v[(v >= tr['t_rs']) & (v < tr['t_rs'] + 0.2)] - tr['t_rs']
            cnt['RS'][k] += np.histogram(q, bins=2, range=(0, 0.2))[0]
            alive = int(np.floor((tr['t_burst'] - tr['t_rs']) / 0.1))
            expo['RS'][k, :min(alive + 1, 2)] += 1
            r = v[(v >= tr['t_resp'] - 0.2) & (v < tr['t_resp'])] - tr['t_resp'] + 0.2
            cnt['response'][k] += np.histogram(r, bins=2, range=(0, 0.2))[0]
            aliveR = int(np.floor((tr['t_burst'] - (tr['t_resp'] - 0.2)) / 0.1))
            expo['response'][k, :max(min(aliveR + 1, 2), 0)] += 1
        for p in PERIODS:
            for k in range(6):
                for b in range(cnt[p].shape[1]):
                    rows.append(dict(subj=s, cell=k, period=p, bin=b,
                                     count=int(cnt[p][k, b]),
                                     ntrial=int(max(expo[p][k, b], 1))))
    return rows


def write(rows, path, fields):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print('wrote', path, len(rows), 'rows')


if __name__ == '__main__':
    write(empirical_rows(), 'pma_counts.csv',
          ['subj', 'cell', 'period', 'bin', 'count', 'ntrial'])
