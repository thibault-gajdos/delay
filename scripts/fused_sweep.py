"""
fused_sweep.py -- the parameter-space evaluation.

What each architecture can produce over its whole admissible parameter space
(manuscript: Methods, "Evaluation over the parameter space"; Results, the "Model--data
comparison" subsections).  Nothing here is fitted.

Latin-hypercube sample over the ranges of Table 1; at every parameter vector, the
admissibility screen and the simulated measures.  The same sampled vector is evaluated in both
models, with beta set to zero in Model 1, so the two are compared at identical parameters.

Ranges (Table 1; their justification is in Methods, "Evaluation over the parameter
space"):

  g        in [0.050, 0.200]
  theta/g  in [0.25, 0.90]                    sampled as a ratio, so that theta < g
  kappa    in [0.20, 30.0]    log-uniform     kappa_epsilon in the manuscript
  SDe      in [0.003, 0.060]  log-uniform     stationary SD of the maintenance noise
  beta     in [0.0015, 0.100] log-uniform     Model 2 only
  k        in [0.5, 2.0]      log-uniform
  lam      in [10.0, 50.0]    log-uniform
  r        in [0.010, 0.200]  log-uniform
  muRS     in [0.30, 1.50]
  sigma_q  in [0.03, 0.30]    log-uniform

Fixed: sigx = 0.100 (sets the amplitude unit), Lmotion = LRS = 0.200 s, delay = 3 s,
Lmotor = 0.093 s.

The cube has 10 dimensions: the nine parameters above shared by both models, and beta,
used in Model 2 only.

Admissibility screen (see admissible()): in every condition, <=10% of trials whose
response burst begins before the response signal, <=10% responses within 0.2 s of the
signal, <=5% no response by the 4 s deadline; and, across conditions, a mean delay-window
partial-activation rate in [0.02, 5] /s.  No bound is derived from the data.

The production run is 1200 draws at 40000 trials per condition per model, simulated in
batches of BATCH trials.  Resumable by draw: a draw already recorded in OUT, or in any file
matching the glob in the environment variable DONE, is skipped.  SHARD_N/SHARD_I in the
environment make shard j take the draws i with i % n == j, each writing its own file; IDX
restricts the run to a JSON list of draws.  The cube is regenerated identically from the
fixed seed in every process, so the record of a draw does not depend on the sharding.

The optional fourth and fifth arguments, GAMMA_X and GAMMA_Q, are the leaks of the
motion-evidence accumulator x and of the response-signal accumulator q (sim.py's
`gamma_x`, `gamma_q`), for the leak-robustness analysis (run_leak.sh).

Usage: python3 fused_sweep.py N_DRAWS NTR OUT_FILE [GAMMA_X [GAMMA_Q]]
"""
import numpy as np, sys, time, json, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS
from repair_shards import repair
from shardio import load

N_DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 250
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = sys.argv[3] if len(sys.argv) > 3 else 'results/fused_sweep_results.jsonl'
GAMMA_X = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
GAMMA_Q = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
SEED0 = 20260824

DELAY_TB = np.arange(20) * 0.1 + 0.05
DELAY_TB = DELAY_TB - DELAY_TB.mean()

RATE_LO, RATE_HI = 0.02, 5.0   # delay-window PMA rate screen, /s (see admissible())

FIXED = dict(sigx=0.100, Lmotion=0.200, LRS=0.200, delay=3.0, Lmotor=0.093)

DIMS = ['g', 'theta_ratio', 'log_kappa', 'log_SDe', 'muRS', 'log_r',
        'log_k', 'log_lam', 'log_sigq', 'log_beta']
LO = dict(g=0.050, theta_ratio=0.25, log_kappa=np.log(0.20), log_SDe=np.log(0.003),
          muRS=0.30, log_r=np.log(0.010),
          log_k=np.log(0.5), log_lam=np.log(10.0), log_sigq=np.log(0.03),
          log_beta=np.log(0.0015))
HI = dict(g=0.200, theta_ratio=0.90, log_kappa=np.log(30.0), log_SDe=np.log(0.060),
          muRS=1.50, log_r=np.log(0.200),
          log_k=np.log(2.0), log_lam=np.log(50.0), log_sigq=np.log(0.30),
          log_beta=np.log(0.100))


def lhs(n, dims, rng):
    """Latin-hypercube sample of n points in len(dims) unit-cube dimensions."""
    d = len(dims)
    cut = np.linspace(0, 1, n + 1)
    u = rng.uniform(size=(n, d))
    samples = np.empty((n, d))
    for j in range(d):
        perm = rng.permutation(n)
        samples[:, j] = cut[perm] + u[:, j] * (cut[perm + 1] - cut[perm])
    return samples


def to_param(row):
    v = dict(zip(DIMS, row))
    g = LO['g'] + v['g'] * (HI['g'] - LO['g'])
    theta_ratio = LO['theta_ratio'] + v['theta_ratio'] * (HI['theta_ratio'] - LO['theta_ratio'])
    kappa = np.exp(LO['log_kappa'] + v['log_kappa'] * (HI['log_kappa'] - LO['log_kappa']))
    SDe = np.exp(LO['log_SDe'] + v['log_SDe'] * (HI['log_SDe'] - LO['log_SDe']))
    muRS = LO['muRS'] + v['muRS'] * (HI['muRS'] - LO['muRS'])
    r = np.exp(LO['log_r'] + v['log_r'] * (HI['log_r'] - LO['log_r']))
    k = np.exp(LO['log_k'] + v['log_k'] * (HI['log_k'] - LO['log_k']))
    lam = np.exp(LO['log_lam'] + v['log_lam'] * (HI['log_lam'] - LO['log_lam']))
    sigq = np.exp(LO['log_sigq'] + v['log_sigq'] * (HI['log_sigq'] - LO['log_sigq']))
    beta_frac = v['log_beta']  # in [0,1]; mapped to alpha below per-model
    theta = g * theta_ratio
    p = dict(FIXED)
    p.update(g=g, theta=theta, kappa=kappa, SDe=SDe, muRS=muRS, r=r,
             k=k, lam=lam, sigq=sigq)
    return p, beta_frac


def counts_profile(mm, t0, t1, nbin):
    """Event counts per bin over [t0, t1) from stimulus onset, and the number of valid
    trials behind them.  The response-related burst is not a PMA and is never counted
    (sim.measures), wherever it begins."""
    v = np.flatnonzero(mm['valid'])
    e = np.concatenate([mm['events'][j] for j in v]) if v.size else np.empty(0)
    h = np.histogram(e[(e >= t0) & (e < t1)], bins=nbin, range=(t0, t1))[0]
    return h.astype(float), v.size


def bin_summary(counts, expo, binw):
    """(level, slope, se_level, se_slope) of a binned rate: the mean rate over the window
    and the least-squares slope against time, in s^-1 and s^-2, each with its Poisson
    Monte-Carlo standard error.  `expo` is the number of trials behind every bin."""
    counts = np.asarray(counts, float)
    expo = np.maximum(np.asarray(expo, float) * np.ones_like(counts), 1.0)
    denom = expo * binw
    r = counts / denom
    n = len(counts)
    x = (np.arange(n) + 0.5) * binw
    x = x - x.mean()
    w1 = x / np.dot(x, x)
    # Poisson: Var(count_i) = count_i, floored at one event, the same convention the
    # delay-window contrasts use
    v = np.maximum(counts, 1.0) / denom ** 2
    return (float(r.mean()), float(w1 @ r),
            float(np.sqrt(np.sum(v)) / n), float(np.sqrt(np.dot(w1 ** 2, v))))


def trend_contrasts(counts, exposure, x):
    """Slope and curvature of the rate against time -- the linear and quadratic
    coefficients of the least-squares fit on the centred bin times in seconds, as in the
    data analysis -- with Poisson Monte-Carlo standard errors.  `counts` are per-bin event
    counts, `exposure` the number of trials backing each bin, `x` the centred bin times.
    Returns (slope, se_slope, curvature, se_curvature) in s^-2 and s^-3."""
    w1 = x / np.dot(x, x)                                   # OLS slope weights
    q = x ** 2 - (x ** 2).mean()
    q = q - x * (np.dot(q, x) / np.dot(x, x))               # orthogonalised
    w2 = q / np.dot(q, q)
    binw = float(np.mean(np.diff(x))) if len(x) > 1 else 1.0
    denom = max(exposure, 1) * binw                         # counts -> rate per second
    r = counts / denom
    lin, quad = float(w1 @ r), float(w2 @ r)
    # Poisson: Var(count_i) = count_i, floored at one event so an empty bin is not
    # credited with zero uncertainty (same convention as the other count measures).
    v = np.maximum(counts, 1.0) / denom ** 2
    return lin, float(np.sqrt(np.dot(w1 ** 2, v))), quad, float(np.sqrt(np.dot(w2 ** 2, v)))


BATCH = 10000   # trials simulated at once; bounds the memory of one process


def simulate_measures(q, c, D, n, model, seed, gamma_x, gamma_q):
    """measures() of n simulated trials of one condition, simulated in batches of at most
    BATCH trials so that the trajectory matrix stays small.  Batch b is seeded with
    (seed, b); the trials of different batches are independent, and the measures are
    concatenated, so the result is a sample of n trials whatever the batch size."""
    parts = []
    for b, start in enumerate(range(0, n, BATCH)):
        nb = min(BATCH, n - start)
        s = simulate(q, c, D, nb, model=model, seed=[seed, b],
                     gamma_x=gamma_x, gamma_q=gamma_q)
        parts.append(measures(s))
        del s
    out = {k: np.concatenate([m[k] for m in parts])
           for k in ('acc', 'rt', 'pmet', 'met', 'valid', 'premature', 'noresp')}
    out['events'] = [e for m in parts for e in m['events']]
    out['D'] = D
    return out


def run_model(p, alpha, model, seed0, gamma_x=0.0, gamma_q=0.0):
    """Simulate all six cells for one (parameter, model) pair.  Returns per-cell arrays
    and the admissibility screen's inputs.  `gamma_x` and `gamma_q` are the leaks of the
    two accumulators (sim.py's `simulate()` arguments of the same names), zero by
    default."""
    q = dict(p)
    q['alpha'] = alpha
    prem = np.zeros(6); early = np.zeros(6); noresp = np.zeros(6); burst_early = np.zeros(6)
    acc = np.zeros(6); rt_mean = np.zeros(6); rt_se = np.zeros(6); rt_n = np.zeros(6)
    rate = np.zeros(6); srate = np.zeros(6)
    trend = np.zeros(6); trend_se = np.zeros(6)
    curv = np.zeros(6); curv_se = np.zeros(6)
    stim_trend = np.zeros(6); stim_trend_se = np.zeros(6)
    met_mean = np.zeros(6); met_se = np.zeros(6)
    pmet_mean = np.zeros(6); pmet_se = np.zeros(6)
    # the 20-bin delay profile per cell, for fig_delay_profiles.py
    prof = np.zeros((6, 20)); prof_expo = np.zeros(6)
    for k, (c, D) in enumerate(CELLS):
        m = simulate_measures(q, c, D, NTR, model, seed0 + k, gamma_x, gamma_q)
        premature = m['premature']
        rt = m['rt']
        early_resp = (~premature) & (~m['noresp']) & (rt < 0.2)
        no_resp = m['noresp'] | ((~premature) & (rt > 4.0))
        prem[k] = premature.mean(); early[k] = early_resp.mean(); noresp[k] = no_resp.mean()
        # trials, out of all of them, whose response burst begins before the response
        # signal (a negative premotor time); every premature response is one
        burst_early[k] = np.mean(np.nan_to_num(m['pmet'], nan=0.0) < 0)
        v = m['valid']
        acc[k] = np.nanmean(m['acc'][v]) if v.sum() else np.nan
        rv = rt[v]
        rt_mean[k] = np.nanmean(rv) if v.sum() else np.nan
        rt_se[k] = np.nanstd(rv) / np.sqrt(max(v.sum(), 1))
        rt_n[k] = v.sum()
        mv = m['met'][v]
        met_mean[k] = np.nanmean(mv) if v.sum() else np.nan
        met_se[k] = np.nanstd(mv) / np.sqrt(max(np.sum(~np.isnan(mv)), 1))
        pv = m['pmet'][v]
        pmet_mean[k] = np.nanmean(pv) if v.sum() else np.nan
        pmet_se[k] = np.nanstd(pv) / np.sqrt(max(np.sum(~np.isnan(pv)), 1))
        # stimulus-onset window: 0 to 2 s from motion onset
        hs, es = counts_profile(m, 0.0, 2.0, 20)
        srate[k] = hs.sum() / (es * 2.0) if es > 0 else np.nan
        _, stim_trend[k], _, stim_trend_se[k] = bin_summary(hs, float(es), 0.1)
        # delay window: 1 to 3 s from motion offset
        h, expo = counts_profile(m, D + 1, D + 3, 20)
        prof[k] = h; prof_expo[k] = expo
        rate[k] = h.sum() / (expo * 2.0) if expo > 0 else np.nan
        if (h > 0).sum() >= 4 and expo > 0:
            trend[k], trend_se[k], curv[k], curv_se[k] = trend_contrasts(h, expo, DELAY_TB)
        else:
            trend[k] = trend_se[k] = curv[k] = curv_se[k] = np.nan
    return dict(prem=prem, early=early, noresp=noresp, burst_early=burst_early,
                acc=acc, rt_mean=rt_mean, rt_se=rt_se, rt_n=rt_n,
                trend=trend, trend_se=trend_se, curv=curv, curv_se=curv_se,
                rate=rate, srate=srate, stim_trend=stim_trend, stim_trend_se=stim_trend_se,
                met_mean=met_mean, met_se=met_se, pmet_mean=pmet_mean, pmet_se=pmet_se,
                prof=prof.astype(np.int64), prof_expo=prof_expo.astype(np.int64))


def admissible(res):
    """A draw stands for a possible participant.  None of the bounds is derived from the
    data.

    In every condition: at most 10% of trials whose response burst begins before the
    response signal, at most 10% responses within 0.2 s of the signal (the response-time
    exclusion rule of the analysis), at most 5% failures to respond by the 4-s deadline
    (the task's own).  The first criterion contains the premature responses, whose burst
    necessarily begins before the signal, and also excludes the regime in which temporal
    preparation alone carries the signal above the gate during the delay, giving one
    sustained activation of several seconds that ends in the response.

    Across conditions: a mean delay-window partial-activation rate in [0.02, 5] /s, which
    excludes regimes with no muscle activity during the delay."""
    if (res['burst_early'] > 0.10).any(): return False
    if (res['early'] > 0.10).any(): return False
    if (res['noresp'] > 0.05).any(): return False
    rm = np.nanmean(res['rate'])
    if not (RATE_LO <= rm <= RATE_HI): return False
    return True


def main():
    rng = np.random.default_rng(SEED0)
    cube = lhs(N_DRAWS, DIMS, rng)
    # Optional sharding. The draws are independent and the LHS cube is
    # regenerated identically from the fixed seed in every process, so a run can be split
    # across cores by having shard j of n take the draws i with i % n == j, each writing
    # its own OUT file; `cat` the shard files to recover the full dump (records carry
    # `draw`, so order is recoverable). Set SHARD_N/SHARD_I in the environment. Default
    # 1/0 = no sharding, byte-identical to the unsharded run for the same draw indices.
    shard_n = int(os.environ.get('SHARD_N', '1'))
    shard_i = int(os.environ.get('SHARD_I', '0'))
    # IDX, if set, is a JSON list of draw indices to evaluate instead of all of them
    # (the leak reruns use the draws admissible at zero leak, from leak_index.py)
    draws = json.load(open(os.environ['IDX'])) if 'IDX' in os.environ else range(N_DRAWS)
    idx = [d for j, d in enumerate(draws) if j % shard_n == shard_i]
    # Resume by draw, not by line count: a draw is skipped if it is already recorded in
    # OUT or in any file matching DONE (a glob, default OUT itself).  The record of a draw
    # depends only on its index, never on the sharding, so the shard count can change
    # between runs without losing or duplicating anything.
    repair(OUT)
    done = set()
    for pat in (os.environ.get('DONE', OUT), OUT):
        for r in load(pat):
            if r.get('ntr') != NTR:
                raise SystemExit(f"{pat} holds draws at another trial count "
                                 f"({r.get('ntr')}, not {NTR}): move it aside first")
            done.add(r['draw'])
    todo = [d for d in idx if d not in done]
    print(f"shard {shard_i}/{shard_n}: {len(idx) - len(todo)} of {len(idx)} draws already "
          f"done (N_DRAWS={N_DRAWS}, NTR={NTR}, gamma_x={GAMMA_X}, gamma_q={GAMMA_Q})",
          flush=True)
    with open(OUT, 'a') as f:
        for i in todo:
            t0 = time.time()
            rec = evaluate(cube, i, GAMMA_X, GAMMA_Q)
            f.write(json.dumps(rec) + "\n"); f.flush(); os.fsync(f.fileno())
            print(f"draw {i}: adm1={rec['model1']['admissible']} "
                  f"adm2={rec['model2']['admissible']} ({time.time()-t0:.1f}s)", flush=True)


def evaluate(cube, i, gamma_x=0.0, gamma_q=0.0):
    """The record of draw i: its parameters, and each model's admissibility and per-cell
    summaries at the given leaks."""
    p, beta_frac = to_param(cube[i])
    beta2 = np.exp(LO['log_beta'] + beta_frac * (HI['log_beta'] - LO['log_beta']))
    # alpha := beta/Delta with Delta = t_RS - t_off = (delay+LRS) - Lmotion,
    # constant across cells since delay/LRS/Lmotion are fixed.
    Delta = (FIXED['delay'] + FIXED['LRS']) - FIXED['Lmotion']
    alpha2 = beta2 / Delta
    res1 = run_model(p, 0.0, 1, SEED0 + 1000 * i, gamma_x, gamma_q)
    res2 = run_model(p, alpha2, 2, SEED0 + 1000 * i + 500, gamma_x, gamma_q)
    rec = dict(draw=i, params=p, beta=beta2, alpha=alpha2, gamma_x=gamma_x,
               gamma_q=gamma_q, ntr=NTR)
    for tag, res in (('model1', res1), ('model2', res2)):
        adm = admissible(res)
        rec[tag] = dict(admissible=adm, raw={k: v.tolist() for k, v in res.items()})
    return rec


if __name__ == '__main__':
    main()
