"""calib_fitB.py -- the all-data fit (fit B).

Fit A (calib_behav_q.py) is the manuscript's primary fit and is untouched: behaviour
only, with the delay-period log-rate slope held out.  Fit B is additional and sits at a
different level -- it asks whether one parameter vector reproduces behaviour and the
electromyographic measures together.  It changes no prediction: the predictions are
derived analytically (Appendices A-E) and tested against the data, not against a fit.

ERROR CHRONOMETRY IS PER CELL (2026-08-31), matching fit A, which scores its incorrect-RT
quantiles per cell at the same error rates.  A version pooling errors across cells was
tried first and is now reverted: it was adopted on the strength of a smoke test at
n = 3000, where the rarest cell yields 3 simulated errors, and that does not describe the
real fit at n = 16000.  Measured there, four of the six cells give 2900 to 8000 simulated
errors and no empty bins; only the two 50%-coherence cells break, with one empty bin each
at about 35 G^2 apiece.  That cost is real but it is not a chronometric failure: at
2.0/c50 the model produces 26 errors in 16000 trials, 0.16% against 5.4% observed, so what
the empty bin registers is the error RATE, which the response-type family already scores.
The double-counting is accepted here because fit A has exactly the same structure -- its
incorrect family is 20.8 of Model 2's 22.7 total -- and consistency between the two fits
matters more than avoiding it.

STRUCTURE.  Ten families per cell (42 categories, against fit A's 13).  Every family is a
multinomial scored by the same deviance g2_term(), imported from calib_behav_q so there is
one implementation:

  response type   {premature, correct, incorrect}          3    <- carries accuracy and
                                                                   the premature rate
  PMET | correct        6 quantile bins                    6
  PMET | incorrect      6                                  6
  MET  | correct        6                                  6
  MET  | incorrect      6                                  6
  delay count     {0, 1, 2, >=3} activations per trial      4
  delay shape     onset over 5 bins of 0.4 s, given >=1     5
  stimulus / RS / response  P(>=1 activation)           2 + 2 + 2

Every family except the first is a CONDITIONAL distribution summing to one on its own.
That is deliberate: the response-type proportions appear in exactly one family, so the
accuracy and premature information is not counted again inside the PMET and MET blocks.
For the same reason RT itself is not a family -- RT is the exact sum PMET + MET (verified
on all 20,562 trials carrying a burst), so scoring all three would triple-count.

WHAT THE OBJECTIVE IS.  A product of marginals, i.e. a composite likelihood, not a single
joint multinomial.  Every caveat attaching to fit A therefore attaches here too and for
the same reason: each family enters at a pseudo sample size of one observation per
participant (N_SUBJ = 22), so the scale of G^2 is a weighting choice rather than the
information the data contain.  No p-value, no confidence interval, no AIC or BIC follows
from it.  With ten families per cell instead of one, fit B's G^2 is on a different scale
from fit A's: THE TWO ARE NOT COMPARABLE TO EACH OTHER.  Only Model 1 against Model 2
within the same fit is meaningful.

MEMORY.  sim.simulate stores the gated drive as an (n x nt) float32 array, about 1 GB at
n = 25000, and measures() then loops over trials in Python.  Simulation is therefore run
in chunks and the statistics accumulated, so n is bounded by time and not by memory.
post = 4.3 s (not sim.simulate's 1.8 default) so that no response inside the 4-s deadline
is truncated into a spurious no-response.

Usage: python3 calib_fitB.py {m1|m2}
"""
import numpy as np, json, os, sys, time
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DEFAULTS
from calib_behav_q import g2_term, bin_fracs, N_SUBJ, FLOOR

TGT = json.load(open('results/targets_b.json'))
TGT_A = json.load(open('results/targets_q.json'))          # premature rate lives here
PREM = {c['label']: c['p_premature'] for c in TGT_A['cells']}
N_SHAPE = TGT['n_shape_bins']
DW = tuple(TGT['delay_window'])
POST = 4.3
CHUNK = 2500

FAMILIES = ['response_type', 'pmet_correct', 'pmet_incorrect', 'met_correct',
            'met_incorrect', 'delay_count', 'delay_shape', 'any_stim', 'any_rs',
            'any_resp']


def sim_stats(p, c, D, n, model, seed, chunk=CHUNK):
    """Model-side counterpart of every fit-B target, accumulated over chunks."""
    q = dict(DEFAULTS); q.update(p)
    delay = q['delay']
    t_rs = D + delay
    pmet_c, pmet_i, met_c, met_i, onsets = [], [], [], [], []
    cnt = []
    any_s = any_r = any_p = 0
    n_valid = n_prem = 0
    done = 0
    while done < n:
        m_ = min(chunk, n - done)
        s = simulate(p, c, D, m_, model=model, seed=seed + done, post=POST)
        mm = measures(s)
        n_prem += int(mm['premature'].sum())
        v = np.flatnonzero(mm['valid'])
        n_valid += v.size
        for j in v:
            ev = mm['events'][j]
            press = t_rs + mm['rt'][j]
            d = ev[(ev >= D + DW[0]) & (ev < D + DW[1])] - (D + DW[0])
            cnt.append(d.size); onsets.extend(d)
            any_s += int(((ev >= 0.0) & (ev < 2.0)).any())
            any_r += int(((ev >= t_rs) & (ev < t_rs + 0.2)).any())
            any_p += int(((ev >= press - 0.2) & (ev < press)).any())
            a, b = mm['pmet'][j], mm['met'][j]
            if np.isfinite(a) and np.isfinite(b):
                if mm['acc'][j]:
                    pmet_c.append(a); met_c.append(b)
                else:
                    pmet_i.append(a); met_i.append(b)
        done += m_
    cnt = np.asarray(cnt, int); onsets = np.asarray(onsets, float)
    nv = max(n_valid, 1)
    if onsets.size:
        h, _ = np.histogram(onsets, bins=N_SHAPE, range=(0.0, DW[1] - DW[0]))
        shape = h / max(h.sum(), 1)
    else:
        shape = np.full(N_SHAPE, 1.0 / N_SHAPE)
    return dict(
        n_prem=n_prem, n_valid=n_valid,
        n_cor=len(pmet_c), n_inc=len(pmet_i),
        pmet_c=np.asarray(pmet_c), pmet_i=np.asarray(pmet_i),
        met_c=np.asarray(met_c), met_i=np.asarray(met_i),
        count_p=np.array([np.mean(cnt == 0), np.mean(cnt == 1),
                          np.mean(cnt == 2), np.mean(cnt >= 3)]) if cnt.size
                 else np.full(4, 0.25),
        shape_p=shape,
        p_any=np.array([any_s / nv, any_r / nv, any_p / nv]),
    )


def cell_g2(cell, st):
    """The ten families for one cell.  Individual family subtotals are not guaranteed
    non-negative -- Gibbs' inequality bounds a complete distribution, and each family is
    complete only within itself, so a subtotal can come out slightly negative when the
    simulated distribution happens to sit closer to uniform than the observed one.  Only
    the sum is the meaningful number."""
    fam = {}

    # 1. response type: premature / correct / incorrect, over premature + valid
    p_prem_o = PREM[cell['label']]
    obs_rt3 = np.array([p_prem_o,
                        (1 - p_prem_o) * cell['p_correct'],
                        (1 - p_prem_o) * cell['p_incorrect']])
    tot = max(st['n_prem'] + st['n_valid'], 1)
    p_prem_s = st['n_prem'] / tot
    nv = max(st['n_valid'], 1)
    sim_rt3 = np.array([p_prem_s,
                        (1 - p_prem_s) * st['n_cor'] / nv,
                        (1 - p_prem_s) * st['n_inc'] / nv])
    fam['response_type'] = (obs_rt3 * N_SUBJ, np.maximum(sim_rt3, FLOOR) * N_SUBJ)

    # 2-5. PMET and MET quantile families, each a conditional distribution
    for name, arr, edges in [
            ('pmet_correct',   st['pmet_c'], cell['pmet_quantiles_correct']),
            ('pmet_incorrect', st['pmet_i'], cell['pmet_quantiles_incorrect']),
            ('met_correct',    st['met_c'],  cell['met_quantiles_correct']),
            ('met_incorrect',  st['met_i'],  cell['met_quantiles_incorrect'])]:
        if not np.all(np.isfinite(edges)):
            continue
        fam[name] = (BIN6 * N_SUBJ, bin_fracs(arr, np.asarray(edges)) * N_SUBJ)

    # 6-7. delay-window activation count and onset shape
    fam['delay_count'] = (np.asarray(cell['delay_count_p']) * N_SUBJ,
                          np.maximum(st['count_p'], FLOOR) * N_SUBJ)
    fam['delay_shape'] = (np.asarray(cell['delay_shape_p']) * N_SUBJ,
                          np.maximum(st['shape_p'], FLOOR) * N_SUBJ)

    # 8-10. P(>=1 activation) in the remaining windows, as two-category families
    for name, key, sp in [('any_stim', 'p_any_stim', st['p_any'][0]),
                          ('any_rs', 'p_any_rs', st['p_any'][1]),
                          ('any_resp', 'p_any_resp', st['p_any'][2])]:
        po = cell[key]
        fam[name] = (np.array([1 - po, po]) * N_SUBJ,
                     np.maximum(np.array([1 - sp, sp]), FLOOR) * N_SUBJ)

    parts = {k: g2_term(o, s) for k, (o, s) in fam.items()}
    return sum(parts.values()), parts


BIN6 = np.array([0.1, 0.2, 0.2, 0.2, 0.2, 0.1])


def objective(p, model, n=10000, detail=False):
    total, parts = 0.0, {}
    for k, cell in enumerate(TGT['cells']):
        st = sim_stats(p, cell['coh'], cell['dur'], n, model, seed=5000 + 100 * k)
        g2, pr = cell_g2(cell, st)
        parts[cell['label']] = pr
        total += g2
    return (total, parts) if detail else total


STAGES = [dict(n=4000, scale=0.10, maxfev=90),
          dict(n=8000, scale=0.05, maxfev=70),
          dict(n=16000, scale=0.02, maxfev=50)]


def run(model, start, tag):
    names = list(start.keys())
    ckpt = f'results/ckptb_{tag}.json'
    stage0 = 0
    if os.path.exists(ckpt):
        st = json.load(open(ckpt)); start = st['p']; stage0 = st.get('stage', 0)
    x0 = np.log([start[k] for k in names])
    t0 = time.time(); p = dict(zip(names, np.exp(x0))); fun = None
    for stage in range(stage0, len(STAGES)):
        cfg = STAGES[stage]
        f = lambda v, n=cfg['n']: objective(dict(zip(names, np.exp(v))), model, n=n)
        sim0 = np.vstack([x0] + [x0 + cfg['scale'] * np.eye(len(x0))[i]
                                 for i in range(len(x0))])
        res = minimize(f, x0, method='Nelder-Mead',
                       options=dict(maxfev=cfg['maxfev'], xatol=1e-4, fatol=1e-3,
                                    adaptive=True, initial_simplex=sim0))
        x0 = res.x; fun = res.fun; p = dict(zip(names, np.exp(x0)))
        json.dump(dict(p=p, obj=float(res.fun), stage=stage + 1, nfev=int(res.nfev)),
                  open(ckpt, 'w'), indent=1)
        with open(f'results/logs/logb_{tag}.txt', 'a') as fh:
            fh.write(f"stage {stage} n={cfg['n']} scale={cfg['scale']}: obj={res.fun:.3f} "
                     f"nfev={res.nfev} elapsed={(time.time()-t0)/60:.1f}min\n")
    return p, fun


if __name__ == '__main__':
    tag = sys.argv[1]
    model = 1 if tag == 'm1' else 2
    start = json.load(open(f'results/calibq_{tag}.json'))   # fit A as the starting point
    if model == 1:
        start.pop('alpha', None)      # inert under Model 1 -- do not search a dead dimension
    p, obj = run(model, start, tag)
    if model == 1:
        p['alpha'] = 0.0
    json.dump(p, open(f'results/fitB_{tag}.json', 'w'), indent=1)
    print(f"{tag}: G2={obj:.3f}")
