"""
export_for_R.py -- the bridge between the Python side and the R side.

Python owns the model, the simulation and the extraction of events from the
electromyographic files; R owns every statistical test.  This script writes the
two files R expects:

    pma_counts.csv      empirical counts, one row per participant, cell, period
                        and bin, with the raw count and the number of trials
                        contributing to that bin.  Needed for the Poisson
                        analyses, which require counts and exposures rather than
                        rates.

    sim_counts_mM.csv   the same quantities from the simulation of model M, with
                        the detection rules applied to the simulated drive, so
                        that the identical R analysis can be run on data and on
                        model output.

Exposure is not always the number of valid trials.  A partial activation is by
definition an activation preceding the response burst, so once that burst has
begun a trial can no longer contribute, and bins after it are undefined rather
than empty.  This matters only in the response-signal and response windows,
where a few per cent of trials have already produced their burst; it is handled
by counting, for each bin, the trials whose burst had not yet begun.

The simulated file carries one participant: a parameter vector is one hypothetical
participant, and the model has no between-participant variability.  The file is
therefore the model's per-cell, per-bin counts and exposures, and supports the
descriptive comparisons the figures draw, not a repeated-measures analysis.
"""
import numpy as np, json, sys, csv, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.empirical import trials
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS

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


def sim_rows(P, model, n_per_cell=22 * 60, seed0=4000):
    rows = []
    for k, (c, D) in enumerate(CELLS):
        s = simulate(P, c, D, n_per_cell, model, seed=seed0 + k)
        m = measures(s)
        # One parameter vector is one simulated participant. An earlier version split a
        # vector's trials into 22 equal groups so that a repeated-measures ANOVA could be
        # run on the model side, but those groups were identical in parameters: their
        # between-group variance was pure sampling noise, and the F values it produced were
        # not comparable with the empirical ones. Nothing in the manuscript reported them.
        valid = np.flatnonzero(m['valid'])
        for g, idx in enumerate([valid]):
            cnt = {p: np.zeros(20 if p in ('stimulus', 'delay') else 2) for p in PERIODS}
            expo = {p: np.zeros(20 if p in ('stimulus', 'delay') else 2) for p in PERIODS}
            for j in idx:
                e = m['events'][j]
                burst = m['pmet'][j] + D + 3.0            # burst onset, from onset
                resp = burst + m['met'][j]
                cnt['stimulus'] += np.histogram(e[(e >= 0) & (e < 2.0)],
                                                bins=20, range=(0, 2))[0]
                expo['stimulus'] += 1
                w = e[(e >= D + 1) & (e < D + 3)] - D - 1
                cnt['delay'] += np.histogram(w, bins=20, range=(0, 2))[0]
                expo['delay'] += 1
                q = e[(e >= D + 3) & (e < D + 3.2)] - D - 3
                cnt['RS'] += np.histogram(q, bins=2, range=(0, 0.2))[0]
                alive = int(np.floor((burst - (D + 3)) / 0.1))
                expo['RS'][:min(alive + 1, 2)] += 1
                r = e[(e >= resp - 0.2) & (e < resp)] - resp + 0.2
                cnt['response'] += np.histogram(r, bins=2, range=(0, 0.2))[0]
                aliveR = int(np.floor((burst - (resp - 0.2)) / 0.1))
                expo['response'][:max(min(aliveR + 1, 2), 0)] += 1
            for p in PERIODS:
                for b in range(len(cnt[p])):
                    nt = max(expo[p][b], 1)
                    rows.append(dict(subj=g, cell=k, period=p, bin=b,
                                     count=int(cnt[p][b]), ntrial=int(nt),
                                     rate=cnt[p][b] / nt))
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
    for mdl in (1, 2):
        R = json.load(open('results/representative.json'))
        P = dict(R['params'], alpha=(R['alpha'] if mdl == 2 else 0.0))
        write(sim_rows(P, mdl), f'sim_counts_m{mdl}.csv',
              ['subj', 'cell', 'period', 'bin', 'count', 'ntrial', 'rate'])
