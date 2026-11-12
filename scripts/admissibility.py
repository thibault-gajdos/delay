"""admissibility.py -- the admissibility screen of the parameter-space evaluation
(manuscript: Methods, "Evaluation over the parameter space").

A parameter vector stands for a possible participant.  It is admissible in a model when,
in every coherence x duration condition (the six cells of targets.CELLS):
  - at least 80% of the trials are retained by the trial exclusions of the data analysis
    (a response before the response signal, no response before the 4-s deadline, or a
    response time below 0.2 s exclude a trial; sim.measures, 'valid');
  - the mean partial-activation rate over the delay window (1 to 3 s after stimulus
    offset, retained trials) lies in [0.02, 5] /s.
Model 2 is moreover restricted to beta < g - theta: in the absence of fluctuations,
temporal preparation does not carry a trial maintained at the internal threshold to the
gate before the response signal.

The screen is applied to the parameters and per-condition quantities stored with every
simulated record (fused_sweep.evaluate), so it can be changed and re-applied without
simulating again.
"""
import numpy as np

MIN_RETAINED = 0.80            # fraction of the trials of a condition retained
RATE_LO, RATE_HI = 0.02, 5.0   # mean delay-window partial-activation rate, /s


def admissible(res, ntr):
    """True when both criteria hold in every condition.  `res` holds per-condition arrays
    (six values each) under the keys 'rt_n' (retained trials) and 'rate'; `ntr` is the
    number of trials simulated per condition."""
    retained = np.asarray(res['rt_n'], float) / ntr
    rate = np.asarray(res['rate'], float)
    if not np.all(retained >= MIN_RETAINED): return False
    if not np.all((rate >= RATE_LO) & (rate <= RATE_HI)): return False   # NaN fails
    return True


def in_space(rec, model):
    """The constraint on the parameter space: beta < g - theta in Model 2; none in Model 1
    (beta = 0)."""
    if model == 'model1':
        return True
    p = rec['params']
    return rec['beta'] < p['g'] - p['theta']


def reapply(rec):
    """Set rec[model]['admissible'] for both models of a simulated record, from its
    parameters and stored per-condition quantities; returns rec."""
    for m in ('model1', 'model2'):
        rec[m]['admissible'] = bool(in_space(rec, m) and admissible(rec[m]['raw'], rec['ntr']))
    return rec
