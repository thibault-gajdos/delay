"""
fused_sweep.py -- the parameter-space evaluation.

What each architecture can produce over its whole admissible parameter space
(manuscript: Methods, "Evaluation over the parameter space"; Results, the "Simulation
validation" subsections).  Nothing here is fitted, so this is where every effect-size
statement in the paper comes from.

Latin-hypercube sample over the ranges of Table 1; at every draw, the task-admissibility
screen and the criteria the Results report.  The same sampled vector is evaluated in both
models, with beta set to zero in Model 1, so the two are compared at identical parameters.

Ranges, and where each comes from.  No bound is derived from these data: a range drawn
around a value estimated from the data cannot test whether a prediction survives away
from that value.

  g        in [0.050, 0.200]                  design: the distance from the gate at theta
                                              to the response bound
  theta/g  in [0.25, 0.90]                    sampled as a ratio, which keeps theta < g
                                              as the model requires
  kappa    in [0.20, 30.0]    log-uniform     1/kappa is the correlation time of the
                                              maintenance noise: [0.20, 30] spans one
                                              trial (5 s) down to one analysis bin
                                              (33 ms), the natural span for a
                                              within-trial fluctuation.  The width
                                              favours neither model: with eps started at
                                              its stationary law the moments of H relax
                                              at rate at least lambda whatever kappa is,
                                              so Model 1's delay rate is flat at every
                                              kappa (Appendix, Corollary A.3).
  SDe      in [0.003, 0.060]  log-uniform     the stationary SD of the maintenance
                                              fluctuation, sampled directly rather than
                                              through the derived a = (g - theta)/SDe:
                                              reparameterising in terms of the quantity
                                              that governs the outcome under test would
                                              redefine the sampling law by its effect on
                                              the result.  The problem that would have
                                              addressed -- draws at a > 3, where
                                              exp(-a^2/2) leaves no crossings at all --
                                              belongs to the acceptance set, and is
                                              handled by the EMG criterion in
                                              admissible().
  beta     in [0.0015, 0.100] log-uniform     Model 2 only.  The upper end is not
                                              derived: the screen enforces it, since a
                                              ramp large enough to carry the signal past
                                              the response criterion during the delay
                                              shows up as premature responses.
  k        in [0.5, 2.0]      log-uniform     about a factor of two either side of
                                              Dendauw et al. (2024)'s 1.183, estimated in
                                              an immediate-report task
  lam      in [10.0, 50.0]    log-uniform     1/lam is the motor-preparation filter time
                                              constant, so the range is tau in [100, 20]
                                              ms; lam = 100 would be 10 ms, faster than
                                              plausible neuromuscular filtering.  This
                                              bound does argumentative work -- lam is
                                              what puts the analysis window in the
                                              stationary regime (Appendix, Proposition
                                              A.1) -- and is stated as a modelling
                                              commitment in the Table 1 note.
  r        in [0.010, 0.200]  log-uniform     response criterion on the neural drive
  muRS     in [0.30, 1.50]                    drift of the response-signal accumulation
  sigma_q  in [0.03, 0.30]    log-uniform     anchored to sigx = 0.100, the other
                                              diffusion SD in the model: sigma_q plays
                                              the analogous role after the response
                                              signal, so a +-3x bracket around sigx is
                                              the least arbitrary anchor available.  The
                                              softest of these bounds.

Fixed: sigx = 0.100 (the model is invariant under a joint rescaling of every amplitude,
so one of them sets the unit), Lmotion = LRS = 0.200 s, delay = 3 s (the task's own), and
Lmotor = 0.093 s (Dendauw et al. 2024 -- a peripheral electromechanical delay has better
reason to transfer between tasks than a decision-related constant, which is why r stays
free).

Latin-hypercube over the 10 dimensions shared by both models (g, theta/g, kappa, SDe,
muRS, r, k, lam, sigma_q and Lmotor's slot); beta is an 11th, sampled for Model 2 only.

Task-admissibility screen (see admissible()): in every condition, <=10% premature, <=10%
responses within 0.2 s of the signal, <=5% no response by the 4 s deadline, accuracy
> 80% in both 50%-coherence conditions; and, across conditions, a delay-window
partial-activation rate in [0.04, 3.0] /s and mean premotor and motor execution times in
the observed per-participant ranges widened by a factor of two.

The production run is 1200 draws at 4000 trials per cell per model per draw (~6.5 s per
draw, ~25 min over six shards).  Resumable by line count of OUT; SHARD_N/SHARD_I in the
environment make shard j take the draws i with i % n == j, each writing its own file.
The cube is regenerated identically from the fixed seed in every process, so cat-ing the
shard files reproduces the unsharded run exactly.

The optional fourth and fifth arguments, GAMMA_X and GAMMA_Q, are the leaks of the
motion-evidence accumulator x and of the response-signal accumulator q (sim.py's
`gamma_x`, `gamma_q`), for the leak-robustness check of the response letter (C2.3);
that check is driven by leak_sweep.py, which reuses this file's cube and screen.

Usage: python3 fused_sweep.py N_DRAWS NTR OUT_FILE [GAMMA_X [GAMMA_Q]]
"""
import numpy as np, sys, time, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS, LAB

N_DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 250
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = sys.argv[3] if len(sys.argv) > 3 else 'results/fused_sweep_results.jsonl'
GAMMA_X = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
GAMMA_Q = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
SEED0 = 20260824

DELAY_TB = np.arange(20) * 0.1 + 0.05
DELAY_TB = DELAY_TB - DELAY_TB.mean()

RATE_LO, RATE_HI = 0.04, 3.0   # delay-window PMA rate screen, /s (see admissible())
# Chronometric screen. See admissible().
PMET_LO, PMET_HI = 0.181, 0.677     # s
MET_LO, MET_HI = 0.046, 0.288       # s

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
    """Raw event counts per bin (not divided by trial count), plus the number of
    valid trials contributing exposure -- constant across delay-window bins, since a
    valid trial's activity is fully observed there (see sim.py measures()).

    One known bias, measured and deliberately left in place.  measures()
    excludes the response-related burst from `events` wherever it begins, because that
    burst IS the response and its onset is what PMET is measured from.  In the model that
    burst sometimes begins BEFORE the response signal -- on up to 2% of trials at some
    admissible settings -- so it is removed from the last two bins of the delay window and
    depresses the window's linear trend.  Participants almost never do this: the empirical
    premotor time is negative on 0.02% of trials.  The bias is about -0.01 s^-2 and applies
    to both architectures, so it makes Model 2's trend conservative and accounts entirely
    for Model 1's simulated median of -0.007 s^-2 against an analytic value of zero.  It is
    left alone rather than patched because the alternative -- counting a burst that became
    the response as a delay-period partial activation -- would be the larger distortion,
    and the correction is three orders of magnitude below the observed trend of +0.56."""
    v = np.flatnonzero(mm['valid'])
    e = np.concatenate([mm['events'][j] for j in v]) if v.size else np.empty(0)
    h = np.histogram(e[(e >= t0) & (e < t1)], bins=nbin, range=(t0, t1))[0]
    return h.astype(float), v.size


def window_stats(mm, D, delay, dt):
    """Counts and per-bin exposures for the three windows outside the delay, following
    exactly the definitions the empirical export uses (export_for_R.empirical_rows):

      stimulus   0 to 2 s from motion onset, 20 bins, every valid trial contributing
      RS         0 to 0.2 s from response-signal onset, 2 bins
      response   0.2 s to 0 s before the button press, 2 bins

    In the last two a trial stops contributing once its response burst has begun, since a
    partial activation is by definition an activation preceding that burst, so the exposure
    is per bin rather than per window.  Times are measured from stimulus onset: the signal
    is at D + delay, the response burst at D + delay + PMET and the response at
    D + delay + RT."""
    v = np.flatnonzero(mm['valid'])
    t_rs = D + delay
    if v.size == 0:
        return dict(stim=(np.zeros(20), 0.0), rs=(np.zeros(2), np.zeros(2)),
                    resp=(np.zeros(2), np.zeros(2)))
    ev = [mm['events'][j] for j in v]
    e = np.concatenate(ev)
    tid = np.repeat(np.arange(v.size), [x.size for x in ev])
    t_burst = t_rs + mm['pmet'][v]
    t_resp = t_rs + mm['rt'][v]
    hs = np.histogram(e[(e >= 0) & (e < 2.0)], bins=20, range=(0, 2))[0]
    q = e[(e >= t_rs) & (e < t_rs + 0.2)] - t_rs
    hrs = np.histogram(q, bins=2, range=(0, 0.2))[0]
    # a trial contributes to bin b of a window while its response burst has not begun,
    # i.e. when the burst onset falls at or after the start of bin b
    alive = np.floor((t_burst - t_rs) / 0.1)
    ers = np.array([(alive >= b).sum() for b in range(2)], float)
    er = e - (t_resp[tid] - 0.2)
    hrp = np.histogram(er[(er >= 0) & (er < 0.2)], bins=2, range=(0, 0.2))[0]
    aliveR = np.floor((t_burst - (t_resp - 0.2)) / 0.1)
    erp = np.array([(aliveR >= b).sum() for b in range(2)], float)
    return dict(stim=(hs.astype(float), float(v.size)), rs=(hrs.astype(float), ers),
                resp=(hrp.astype(float), erp))


def bin_summary(counts, expo, binw):
    """(level, trend, se_level, se_trend) of a binned rate: the mean rate over the window
    and the linear orthogonal-polynomial contrast against time, in s^-1 and s^-2, each
    with its Poisson Monte-Carlo standard error.  `expo` may be a scalar (constant
    exposure) or one value per bin, which is what the response-signal and response windows
    need, since a trial stops contributing once its response burst has begun."""
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
    """Linear and quadratic orthogonal-polynomial contrasts of the RAW rate against
    time, with Poisson Monte-Carlo standard errors.  This is the statistic the data
    analysis uses (Methods: TREND and QUAD on the binned rate); the offset-Poisson
    log-rate slope below is a different quantity and is kept only as a secondary
    record.  `counts` are per-bin event counts, `exposure` the number of trials
    backing each bin, `x` the centred bin times in seconds.  Returns
    (linear, se_linear, quadratic, se_quadratic) in s^-2 and s^-3."""
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


def poisson_slope_se(counts, exposure, x, iters=100):
    """Offset Poisson-GLM slope (log mu = log(exposure) + b0 + b1*x) and its analytic
    SE from the final Fisher information. `exposure` (the number of trials backing
    each bin) MUST enter as an offset: fitting the already-normalized rate as if it
    were itself a Poisson count (dividing out `exposure` first) discards the sample
    size and makes the resulting SE far too large, independent of how many trials
    were simulated."""
    off = np.log(max(exposure, 1))
    X = np.column_stack([np.ones_like(x), x])
    y = counts
    b = np.array([np.log(max(y.mean(), 1e-4)) - off, 0.0])
    A = None
    for _ in range(iters):
        mu = np.exp(np.clip(off + X @ b, -30, 30))
        z = X @ b + (y - mu) / np.maximum(mu, 1e-12)
        A = X.T @ (mu[:, None] * X) + 1e-9 * np.eye(2)
        bn = np.linalg.solve(A, X.T @ (mu * z))
        if np.max(np.abs(bn - b)) < 1e-11:
            b = bn
            break
        b = bn
    cov = np.linalg.inv(A)
    return b[1], np.sqrt(cov[1, 1])



def run_model(p, alpha, model, seed0, gamma_x=0.0, gamma_q=0.0):
    """Simulate all six cells for one (parameter, model) pair.  Returns per-cell arrays
    and admissibility/criterion summaries.  `gamma_x` and `gamma_q` are the leaks of the
    two accumulators (sim.py's `simulate()` arguments of the same names), zero by
    default."""
    q = dict(p)
    q['alpha'] = alpha
    prem = np.zeros(6); early = np.zeros(6); noresp = np.zeros(6)
    acc = np.zeros(6); rt_mean = np.zeros(6); rt_se = np.zeros(6); rt_n = np.zeros(6)
    slope = np.zeros(6); slope_se = np.zeros(6); rate = np.zeros(6)
    trend = np.zeros(6); trend_se = np.zeros(6)
    curv = np.zeros(6); curv_se = np.zeros(6)
    # the raw 20-bin delay profile per cell, so that figures_space.py can draw the
    # simulated delay-period rate curve against the observed one without re-simulating
    prof = np.zeros((6, 20)); prof_expo = np.zeros(6)
    met_mean = np.zeros(6); met_se = np.zeros(6)
    pmet_mean = np.zeros(6); pmet_se = np.zeros(6)
    # stimulus-period partial-activation rate, 0-2 s from motion onset, the window
    # empirical.py uses; carries Proposition B.2, the only prediction outside the delay
    # on which the two models disagree
    srate = np.zeros(6)
    # level and linear trend in the three windows outside the delay, so that the effects
    # the Results report there can be compared with the model as the delay-period ones are
    stim_trend = np.zeros(6); stim_trend_se = np.zeros(6)
    rs_rate = np.zeros(6); rs_trend = np.zeros(6)
    rs_rate_se = np.zeros(6); rs_trend_se = np.zeros(6)
    resp_rate = np.zeros(6); resp_trend = np.zeros(6)
    resp_rate_se = np.zeros(6); resp_trend_se = np.zeros(6)
    for k, (c, D) in enumerate(CELLS):
        s = simulate(q, c, D, NTR, model=model, seed=seed0 + k,
                     gamma_x=gamma_x, gamma_q=gamma_q)
        m = measures(s)
        ri = s['resp_i']; i_phys = int(round((D + q['delay']) / DT))
        premature = (ri >= 0) & (ri < i_phys)
        rt = m['rt']
        early_resp = (~premature) & (ri >= 0) & (rt < 0.2)
        no_resp = (ri < 0) | ((~premature) & (rt > 4.0))
        prem[k] = premature.mean(); early[k] = early_resp.mean(); noresp[k] = no_resp.mean()
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
        hs, es = counts_profile(m, 0.0, 2.0, 20)
        srate[k] = hs.sum() / (es * 2.0) if es > 0 else np.nan
        # the three windows outside the delay, on the definitions the data analysis uses
        ws = window_stats(m, D, q['delay'], DT)
        _, stim_trend[k], _, stim_trend_se[k] = bin_summary(*ws['stim'], 0.1)
        rs_rate[k], rs_trend[k], rs_rate_se[k], rs_trend_se[k] = bin_summary(*ws['rs'], 0.1)
        resp_rate[k], resp_trend[k], resp_rate_se[k], resp_trend_se[k] = \
            bin_summary(*ws['resp'], 0.1)
        h, expo = counts_profile(m, D + 1, D + 3, 20)
        prof[k] = h; prof_expo[k] = expo
        # mean delay-period partial-activation rate over the 2-s analysis window, per
        # trial per second -- directly comparable with the observed end-of-delay rate.
        rate[k] = h.sum() / (expo * 2.0) if expo > 0 else np.nan
        if (h > 0).sum() >= 4 and expo > 0:
            slope[k], slope_se[k] = poisson_slope_se(h, expo, DELAY_TB)
            trend[k], trend_se[k], curv[k], curv_se[k] = trend_contrasts(h, expo, DELAY_TB)
        else:
            slope[k], slope_se[k] = np.nan, np.nan
            trend[k] = trend_se[k] = curv[k] = curv_se[k] = np.nan
    return dict(prem=prem, early=early, noresp=noresp, acc=acc, rt_mean=rt_mean,
                rt_se=rt_se, rt_n=rt_n, slope=slope, slope_se=slope_se,
                trend=trend, trend_se=trend_se, curv=curv, curv_se=curv_se,
                rate=rate, srate=srate,
                stim_trend=stim_trend, stim_trend_se=stim_trend_se,
                rs_rate=rs_rate, rs_rate_se=rs_rate_se,
                rs_trend=rs_trend, rs_trend_se=rs_trend_se,
                resp_rate=resp_rate, resp_rate_se=resp_rate_se,
                resp_trend=resp_trend, resp_trend_se=resp_trend_se,
                met_mean=met_mean, met_se=met_se,
                pmet_mean=pmet_mean, pmet_se=pmet_se,
                prof=prof.astype(np.int64), prof_expo=prof_expo.astype(np.int64))


def admissible(res):
    if (res['prem'] > 0.10).any(): return False
    if (res['early'] > 0.10).any(): return False
    if (res['noresp'] > 0.05).any(): return False
    # cells index 2 and 5 are c=0.50, D=0.3 and D=2.0 (see targets.CELLS)
    if not (res['acc'][2] > 0.80 and res['acc'][5] > 0.80): return False
    # EMG criterion.  Every criterion above is behavioural, and a purely behavioural
    # screen admits settings that produce literally zero delay-period muscle activity --
    # about half the behaviourally admissible Model 1 draws.  A draw of parameters stands for a
    # hypothetical participant, so the criterion is stated at the participant level: the
    # draw's mean delay-window partial-activation rate must fall inside the range the
    # real participants span.  Observed: 132 subject x condition cells run 0.085-2.785
    # per second (mean 0.709); per-participant means run 0.180-1.863.  The window used
    # here, RATE_LO..RATE_HI = 0.04..3.0, is deliberately WIDER than either -- the
    # conclusion this sweep supports is an EXCLUSION of Model 1, so erring wide is the
    # conservative direction.  In practice the ceiling never binds (no admissible draw of
    # either model exceeded 3.0); the criterion is a floor.
    rm = np.nanmean(res['rate'])
    if not (RATE_LO <= rm <= RATE_HI): return False
    # Chronometric criterion, same form and same justification as the
    # rate window above: a draw stands for a hypothetical participant, and a participant
    # whose response-time decomposition is physically impossible is not one. Without it the
    # admissible set contained draws with a mean premotor execution time of -945 ms -- a
    # response burst beginning almost a second BEFORE the reaction signal -- and mean motor
    # execution times of 1.4 s.
    #
    # The window is the range of the observed per-participant means, widened by a factor of
    # two about their mean: observed PMET means run 270-518 ms and MET means 103-224 ms,
    # giving 181-677 ms and 46-288 ms here.
    #
    # This screens on the LEVELS of the two components, and it is important that it does not
    # screen on their EFFECTS. The paper's RS-locked predictions are that coherence and
    # duration move the premotor component and leave the motor component alone; those are
    # statements about differences between conditions, and they remain free to be right or
    # wrong at every admissible draw. Screening on the level of PMET therefore removes no
    # tested prediction, where screening on the invariance of MET would have converted a
    # result into an assumption -- which is why that was considered and rejected.
    pm = np.nanmean(res['pmet_mean']); mt = np.nanmean(res['met_mean'])
    if not (np.isfinite(pm) and PMET_LO <= pm <= PMET_HI): return False
    if not (np.isfinite(mt) and MET_LO <= mt <= MET_HI): return False
    return True


def rt_effects(res):
    """Coherence slope (marginal over duration, 3 levels) and duration difference
    (marginal over coherence), with MC-SE, mirroring the manuscript's existing
    RT-effect test (para 86)."""
    rt = res['rt_mean']; se = res['rt_se']
    coh_marg = np.array([(rt[0] + rt[3]) / 2, (rt[1] + rt[4]) / 2, (rt[2] + rt[5]) / 2])
    coh_se = np.array([np.sqrt(se[0]**2 + se[3]**2) / 2,
                        np.sqrt(se[1]**2 + se[4]**2) / 2,
                        np.sqrt(se[2]**2 + se[5]**2) / 2])
    x = np.array([0.0, 0.12, 0.50]); x = x - x.mean()
    w = x / (x @ x)
    coh_slope = w @ coh_marg
    coh_slope_se = np.sqrt((w**2) @ (coh_se**2))
    dur_short = rt[0:3].mean(); dur_long = rt[3:6].mean()
    dur_short_se = np.sqrt((se[0:3]**2).sum()) / 3
    dur_long_se = np.sqrt((se[3:6]**2).sum()) / 3
    dur_diff = dur_short - dur_long
    dur_diff_se = np.sqrt(dur_short_se**2 + dur_long_se**2)
    coh_supported = (coh_slope + 1.96 * coh_slope_se) < 0
    dur_supported = (dur_diff - 1.96 * dur_diff_se) > 0
    return dict(coh_slope=coh_slope, coh_slope_se=coh_slope_se,
                dur_diff=dur_diff, dur_diff_se=dur_diff_se,
                coh_supported=bool(coh_supported), dur_supported=bool(dur_supported))


def slope_criterion(res):
    """All six delay-period log-rate slopes significantly positive (lower 95% CI > 0)."""
    lo = res['slope'] - 1.96 * res['slope_se']
    ok = np.isfinite(lo).all() and (lo > 0).all()
    return dict(all_positive=bool(ok), min_slope=float(np.nanmin(res['slope'])),
                max_slope=float(np.nanmax(res['slope'])))


def trend_criterion(res):
    """All six delay-period raw-rate TRENDS significantly positive (lower 95% CI > 0).

    The counterpart of slope_criterion, and the statistic the Results cite: the trend is
    linear in the rate, so unlike a log-rate slope it cannot be earned by ramping from an
    arbitrarily low level.  Both are recorded, and they agree."""
    lo = res['trend'] - 1.96 * res['trend_se']
    ok = bool(np.isfinite(lo).all() and (lo > 0).all())
    return dict(all_positive=ok, min_trend=float(np.nanmin(res['trend'])),
                max_trend=float(np.nanmax(res['trend'])))


def met_invariance(res):
    """Range of per-cell mean MET: the model predicts it invariant across conditions,
    to within the few milliseconds the filter transient leaves."""
    rng = float(np.nanmax(res['met_mean']) - np.nanmin(res['met_mean']))
    return dict(range_ms=rng * 1000.0)


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
    idx = [i for i in range(N_DRAWS) if i % shard_n == shard_i]
    if os.path.exists(OUT):
        done = sum(1 for _ in open(OUT))
    else:
        done = 0
    print(f"resuming at {done}/{len(idx)} of shard {shard_i}/{shard_n} "
          f"(N_DRAWS={N_DRAWS}, NTR={NTR})", flush=True)
    with open(OUT, 'a') as f:
        for i in idx[done:]:
            t0 = time.time()
            p, beta_frac = to_param(cube[i])
            beta2 = np.exp(LO['log_beta'] + beta_frac *
                           (HI['log_beta'] - LO['log_beta']))
            # alpha := beta/Delta with Delta = t_RS - t_off = (delay+LRS) - Lmotion,
            # constant across cells since delay/LRS/Lmotion are fixed.
            Delta = (FIXED['delay'] + FIXED['LRS']) - FIXED['Lmotion']
            alpha2 = beta2 / Delta
            res1 = run_model(p, 0.0, 1, SEED0 + 1000 * i, GAMMA_X, GAMMA_Q)
            res2 = run_model(p, alpha2, 2, SEED0 + 1000 * i + 500, GAMMA_X, GAMMA_Q)
            adm1 = admissible(res1); adm2 = admissible(res2)
            rec = dict(draw=i, params=p, beta=beta2, alpha=alpha2, gamma_x=GAMMA_X, gamma_q=GAMMA_Q,
                       model1=dict(admissible=adm1,
                                   rt=rt_effects(res1) if adm1 else None,
                                   slope=slope_criterion(res1) if adm1 else None,
                                   met=met_invariance(res1) if adm1 else None,
                                   raw={k: v.tolist() for k, v in res1.items()}),
                       model2=dict(admissible=adm2,
                                   rt=rt_effects(res2) if adm2 else None,
                                   slope=slope_criterion(res2) if adm2 else None,
                                   met=met_invariance(res2) if adm2 else None,
                                   raw={k: v.tolist() for k, v in res2.items()}))
            f.write(json.dumps(rec) + "\n"); f.flush()
            print(f"draw {i}: adm1={adm1} adm2={adm2} ({time.time()-t0:.1f}s)", flush=True)


if __name__ == '__main__':
    main()
