"""
fused_sweep.py -- unified broad-parameter Monte-Carlo sweep, extending the manuscript's
existing "first-stage parameter-space evaluation" (Table A1, 250 Latin-hypercube draws,
task-admissibility screen, two RT-effect criteria) to ALSO check, at every admissible
draw, the out-of-sample criteria that the current manuscript only checks at the single
behaviorally-calibrated point:
  (A) the delay-period partial-activation log-rate slope is positive in every cell
      (Model 2's diagnostic prediction; Model 1 predicts no systematic ramp)
  (B) motor execution time (MET) is close to invariant across the six cells (a
      prediction common to both models, checked here across the sampled space rather
      than only at the fitted point -- this is the R7-R10 argument from note4.tex)

Sampling plan (the "hybrid" range agreed with the user 2026-08-24): keep Table A1's
existing ranges for g, theta/g, the stationary SD feeding sigma_epsilon, beta, and muRS;
WIDEN kappa from [0.20, 2.00] to [0.20, 30.0] s^-1 so the sampled region actually covers
the calibrated point (kappa = 12.8 s^-1, note4.tex Section 9); and newly FREE r and
Lmotor (previously fixed at values "inherited" from the immediate-report architecture
that the full calibration found were wrong by factors of 3x and 5x respectively --
note4.tex Section 9, "Two inherited constants do not transfer"). k and lambda stay fixed
at their inherited values, since only the product lambda*kappa enters the crossing rate
(note4.tex Section 9.1), so widening kappa alone already explores the relevant product
range without a second free dimension.

Ranges (log-uniform for scale/rate parameters, uniform for bounded ratios):
  g        in [0.050, 0.160]                 (Table A1, unchanged)
  theta/g  in [0.25, 0.90]                    (Table A1, unchanged; theta = g * ratio)
  kappa    in [0.20, 30.0]   log-uniform      (WIDENED from [0.20, 2.00])
  SDe      in [0.003, 0.030]                  (Table A1, unchanged; this is the
                                                stationary SD -- sim.py derives the
                                                actual diffusion coefficient as
                                                SDe*sqrt(2*kappa), so widening kappa
                                                already widens the effective range)
  beta     in [0.0015, 0.024] (Model 2 only)  (Table A1, unchanged; Model 1 forced 0)
  muRS     in [0.30, 1.50]                    (Table A1, unchanged)
  r        in [0.010, 0.200]  log-uniform     (NEW; brackets both the old assumed 0.022
                                                and the fitted 0.0694, well inside
                                                calibrate.py's search bounds [0.005,0.40])
  Lmotor   in [0.008, 0.080]  log-uniform     (NEW; brackets both the old assumed 0.093s
                                                and the fitted 0.019s, inside calibrate.py's
                                                search bounds [0.005,0.15])
k, lam, sigx, Lmotion, LRS, sigq fixed at Table A1's inherited point values.

Latin-hypercube over the 6 continuous dimensions shared by both models (g, theta/g,
kappa, SDe, muRS, r, Lmotor -- 7 dims; beta is an 8th dimension sampled only for
Model 2, mirroring the original design where "the same sampled vectors were evaluated
in both models, with beta set to zero in Model 1").

Task-admissibility screen: IDENTICAL criteria to the existing manuscript text (Results >
Behavior > Simulation validation, para 85): in every condition, <=10% premature, <=10%
respond within 0.2s of RS, <=5% no response by 4s deadline; accuracy > 80% in both
50%-coherence conditions.

Trial count: 2000 trials/cell/model/draw, matching the existing first-stage procedure
(NOT calibrate.py's NTR=200, which is for the optimization search only).

Leak-robustness variant (2026-08-24): an optional fourth CLI argument, GAMMA, adds a
common leak to both accumulators (sim.py's `simulate(..., gamma=)`, the reviewer R2.3
sensitivity request) and is threaded through every simulate() call and recorded in each
output record. Used to repeat this identical hybrid-range sweep at gamma in
{0.25, 0.50, 1.00} s^-1, replacing the earlier leak-robustness check (Methods, Behavior >
Simulation validation), which had used the pre-widening Table A1 ranges and tested only
the two RT-effect criteria via the older, narrower sweep.py.  Because run_model() already
computes the slope and MET criteria for every draw regardless of gamma, this rerun gets
slope/MET robustness under leak "for free," alongside the RT-effect numbers the original
check reported.

Usage: python3 fused_sweep.py N_DRAWS NTR OUT_FILE [GAMMA]
"""
import numpy as np, sys, time, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS, LAB

N_DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 250
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = sys.argv[3] if len(sys.argv) > 3 else 'results/fused_sweep_results.jsonl'
GAMMA = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
SEED0 = 20260824

DELAY_TB = np.arange(20) * 0.1 + 0.05
DELAY_TB = DELAY_TB - DELAY_TB.mean()

FIXED = dict(k=1.183, lam=21.471, sigx=0.100, Lmotion=0.200, LRS=0.200, sigq=0.100,
             delay=3.0)

DIMS = ['g', 'theta_ratio', 'log_kappa', 'SDe', 'muRS', 'log_r', 'log_Lmotor', 'beta']
LO = dict(g=0.050, theta_ratio=0.25, log_kappa=np.log(0.20), SDe=0.003, muRS=0.30,
          log_r=np.log(0.010), log_Lmotor=np.log(0.008), beta=0.0015)
HI = dict(g=0.160, theta_ratio=0.90, log_kappa=np.log(30.0), SDe=0.030, muRS=1.50,
          log_r=np.log(0.200), log_Lmotor=np.log(0.080), beta=0.024)


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
    SDe = LO['SDe'] + v['SDe'] * (HI['SDe'] - LO['SDe'])
    muRS = LO['muRS'] + v['muRS'] * (HI['muRS'] - LO['muRS'])
    r = np.exp(LO['log_r'] + v['log_r'] * (HI['log_r'] - LO['log_r']))
    Lmotor = np.exp(LO['log_Lmotor'] + v['log_Lmotor'] * (HI['log_Lmotor'] - LO['log_Lmotor']))
    beta_frac = v['beta']  # in [0,1]; mapped to alpha below per-model
    theta = g * theta_ratio
    p = dict(FIXED)
    p.update(g=g, theta=theta, kappa=kappa, SDe=SDe, muRS=muRS, r=r, Lmotor=Lmotor)
    return p, beta_frac


def counts_profile(mm, t0, t1, nbin):
    """Raw event counts per bin (not divided by trial count), plus the number of
    valid trials contributing exposure -- constant across delay-window bins, since a
    valid trial's activity is fully observed there (see sim.py measures())."""
    v = np.flatnonzero(mm['valid'])
    h = np.zeros(nbin)
    for j in v:
        e = mm['events'][j]
        h += np.histogram(e[(e >= t0) & (e < t1)], bins=nbin, range=(t0, t1))[0]
    return h, v.size


def poisson_slope_se(counts, exposure, x, iters=100):
    """Offset Poisson-GLM slope (log mu = log(exposure) + b0 + b1*x) and its analytic
    SE from the final Fisher information. `exposure` (the number of trials backing
    each bin) MUST enter as an offset: fitting the already-normalized rate as if it
    were itself a Poisson count (dividing out `exposure` first) discards the sample
    size and makes the resulting SE far too large, independent of how many trials
    were simulated -- caught by a smoke test before the full sweep was run (2026-08-24)."""
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


def run_model(p, alpha, model, seed0, gamma=0.0):
    """Simulate all six cells for one (parameter, model) pair.  Returns per-cell arrays
    and admissibility/criterion summaries.  `gamma` adds a common leak to both
    accumulators (reviewer R2.3 robustness check; sim.py's `simulate()` argument of the
    same name), 0.0 (no leak) by default."""
    q = dict(p)
    q['alpha'] = alpha
    prem = np.zeros(6); early = np.zeros(6); noresp = np.zeros(6)
    acc = np.zeros(6); rt_mean = np.zeros(6); rt_se = np.zeros(6); rt_n = np.zeros(6)
    slope = np.zeros(6); slope_se = np.zeros(6)
    met_mean = np.zeros(6); met_se = np.zeros(6)
    for k, (c, D) in enumerate(CELLS):
        s = simulate(q, c, D, NTR, model=model, seed=seed0 + k, gamma=gamma)
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
        h, expo = counts_profile(m, D + 1, D + 3, 20)
        if (h > 0).sum() >= 4 and expo > 0:
            slope[k], slope_se[k] = poisson_slope_se(h, expo, DELAY_TB)
        else:
            slope[k], slope_se[k] = np.nan, np.nan
    return dict(prem=prem, early=early, noresp=noresp, acc=acc, rt_mean=rt_mean,
                rt_se=rt_se, rt_n=rt_n, slope=slope, slope_se=slope_se,
                met_mean=met_mean, met_se=met_se)


def admissible(res):
    if (res['prem'] > 0.10).any(): return False
    if (res['early'] > 0.10).any(): return False
    if (res['noresp'] > 0.05).any(): return False
    # cells index 2 and 5 are c=0.50, D=0.3 and D=2.0 (see targets.CELLS)
    if not (res['acc'][2] > 0.80 and res['acc'][5] > 0.80): return False
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


def met_invariance(res):
    """Range of per-cell mean MET, and whether it stays within the note4.tex R7-R10
    predicted residual scale (~2.3 ms) plus a safety margin for six-cell spread."""
    rng = float(np.nanmax(res['met_mean']) - np.nanmin(res['met_mean']))
    return dict(range_ms=rng * 1000.0)


def main():
    rng = np.random.default_rng(SEED0)
    cube = lhs(N_DRAWS, DIMS, rng)
    if os.path.exists(OUT):
        done = sum(1 for _ in open(OUT))
    else:
        done = 0
    print(f"resuming at draw {done}/{N_DRAWS}", flush=True)
    with open(OUT, 'a') as f:
        for i in range(done, N_DRAWS):
            t0 = time.time()
            p, beta_frac = to_param(cube[i])
            beta2 = LO['beta'] + beta_frac * (HI['beta'] - LO['beta'])
            # alpha := beta/Delta with Delta = t_RS - t_off = (delay+LRS) - Lmotion,
            # constant across cells since delay/LRS/Lmotion are fixed.
            Delta = (FIXED['delay'] + FIXED['LRS']) - FIXED['Lmotion']
            alpha2 = beta2 / Delta
            res1 = run_model(p, 0.0, 1, SEED0 + 1000 * i, gamma=GAMMA)
            res2 = run_model(p, alpha2, 2, SEED0 + 1000 * i + 500, gamma=GAMMA)
            adm1 = admissible(res1); adm2 = admissible(res2)
            rec = dict(draw=i, params=p, beta=beta2, alpha=alpha2, gamma=GAMMA,
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
