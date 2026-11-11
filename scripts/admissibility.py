"""admissibility.py -- the admissibility screen of the parameter-space evaluation
(manuscript: Methods, "Evaluation over the parameter space").

A parameter vector stands for a possible participant.  It is admissible in a model when,
in every coherence x duration condition (the six cells of targets.CELLS):
  - a response occurs before the response signal on at most 10% of trials;
  - a response occurs within 0.2 s after the response signal on at most 10% of trials;
  - no response occurs before the 4-s deadline on at most 5% of trials;
  - the mean partial-activation rate over the delay window (1 to 3 s after stimulus
    offset) lies in [0.02, 5] /s.
All fractions are taken over all the trials of the condition.

The screen is applied to the per-condition quantities stored with every simulated record
(fused_sweep.run_model), so it can be changed and re-applied without simulating again.
"""
import numpy as np

MAX_PREMATURE = 0.10      # response before the response signal
MAX_EARLY = 0.10          # response within 0.2 s after the response signal
MAX_NORESP = 0.05         # no response before the 4-s deadline
RATE_LO, RATE_HI = 0.02, 5.0   # mean delay-window partial-activation rate, /s


def admissible(res):
    """True when every criterion holds in every condition.  `res` holds per-condition
    arrays (six values each) under the keys 'prem', 'early', 'noresp' and 'rate'."""
    prem = np.asarray(res['prem'], float)
    early = np.asarray(res['early'], float)
    noresp = np.asarray(res['noresp'], float)
    rate = np.asarray(res['rate'], float)
    if not np.all(prem <= MAX_PREMATURE): return False
    if not np.all(early <= MAX_EARLY): return False
    if not np.all(noresp <= MAX_NORESP): return False
    if not np.all((rate >= RATE_LO) & (rate <= RATE_HI)): return False   # NaN fails
    return True


def reapply(rec):
    """Set rec[model]['admissible'] from the stored per-condition quantities, for both
    models of a simulated record; returns rec."""
    for m in ('model1', 'model2'):
        rec[m]['admissible'] = bool(admissible(rec[m]['raw']))
    return rec
