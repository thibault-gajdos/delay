"""
calibrate.py -- shared simulation/statistics utility, `model_stats()`.

Removed 2026-08-30: the joint EMG+behaviour fit that used to live here (Nelder-Mead
search over all thirteen free parameters against nine blocks per cell -- accuracy,
premotor/motor execution time and their spread, delay-period rate and log-rate slope,
stimulus-period rate, premature-response rate -- each standardized by its
between-participant SE). That fit's own match to the delay-period rate slope was
necessarily in-sample, since `slope` was one of the nine blocks in its objective --
it could not serve as a validation of the model's electromyographic predictions.

The project's calibration is now behaviour alone (`calib_behav_q.py`, an RT-quantile
G² criterion with no EMG term at all -- see its own docstring and the README), and the
electromyographic comparison is instead an honest out-of-sample check: simulate the
behaviour-only-calibrated parameters and compare their predicted delay-period rate
profile against the real EMG (`figures.py`'s `_oos_predictions`/`fig_oos_*`, and
`results/oos_slope_check.json`).

`model_stats()` below is retained because it is the one place that turns a parameter
dict into the full set of per-cell simulated quantities -- accuracy, execution times,
delay/stimulus rate, premature rate, the delay log-rate slope and the full rate
profile -- and both `sweep.py` (robustness table) and the out-of-sample check need it.
It computes; it does not fit.

Recover the removed fitting code (`make_objective`, `ALL_BLOCKS`/`BEHAV_BLOCKS`,
`START`/`LO`/`HI`/`clip`/`pack`/`unpack`, the `__main__` optimizer) from git history if
the question of fitting the model to actual EMG (using the same G²-style methods as
the behaviour-only fit, once we've decided what belongs in that G²) comes back --
`git log --all -p -- calibrate.py` around the 2026-08-30 commits.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, rate_profile
from delaymodel.targets import CELLS, poisson_slope

NTR = 200
TB = np.arange(20) * 0.1 + 0.05; TB = TB - TB.mean()

STAT_KEYS = ['acc', 'pmet', 'met', 'pmet_sd', 'met_sd', 'dly', 'stim', 'prem', 'slope']


def model_stats(p, n=NTR, model=2, seed0=1000, gamma=0.0):
    out = {k: np.zeros(6) for k in STAT_KEYS}
    prof = np.zeros((6, 20))
    out['rt'] = np.zeros(6); out['rt_sd'] = np.zeros(6)
    for k, (c, D) in enumerate(CELLS):
        s = simulate(p, c, D, n, model, seed=seed0 + k, gamma=gamma)
        m = measures(s)
        v = m['valid']
        out['prem'][k] = m['premature'].mean()
        if v.sum() < 10:
            out['acc'][k] = np.nan; continue
        out['acc'][k] = m['acc'][v].mean()
        rv = m['rt'][v]
        out['rt'][k] = np.nanmean(rv)
        out['rt_sd'][k] = np.subtract(*np.nanpercentile(rv, [75, 25]))
        out['pmet'][k] = np.nanmean(m['pmet'][v]); out['met'][k] = np.nanmean(m['met'][v])
        pv = m['pmet'][v]; mv = m['met'][v]
        out['pmet_sd'][k] = np.subtract(*np.nanpercentile(pv, [75, 25]))
        out['met_sd'][k] = np.subtract(*np.nanpercentile(mv, [75, 25]))
        pr = rate_profile(m, D + 1, D + 3, 20); prof[k] = pr
        out['dly'][k] = pr.mean()
        out['stim'][k] = rate_profile(m, 0.0, 2.0, 20).mean()
        out['slope'][k] = poisson_slope(pr, TB)[1] if (pr > 0).sum() >= 4 else np.nan
    out['profile'] = prof
    return out
