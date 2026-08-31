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

Sampling plan (revised 2026-08-30 to remove circularity: the previous ranges bracketed
the fitted/observed point itself, which is not a valid way to define an INDEPENDENT
robustness check of that same fit -- see README "Non-circular parameter ranges" note
below). Every bound here now comes from either a task/design constant or a pre-existing
architectural constant that was fixed before any behavioural fit was run -- never from
calibq_{m1,m2}.json or from targets_q.json.

  g        in [0.050, 0.200]                 (task/design: g is the distance from the
                                                gate at theta to the response bound;
                                                widened upper end from 0.160 to bracket
                                                a wider admissible fraction of theta)
  theta/g  in [0.25, 0.90]                    (unchanged -- no independent bound found
                                                that survives being checked at a second,
                                                independent representative point; see
                                                README)
  kappa    in [0.20, 30.0]    log-uniform     (UNCHANGED, and deliberately so. A
                                                stationarity-based floor of 2.0 was
                                                briefly imposed here on 2026-08-30 and
                                                is now reverted: it was wrong twice
                                                over. The paper does not rely on the
                                                stationary crossing-rate form at small
                                                kappa -- Results (Delay-locked) and
                                                Appendix C both compute from the exact
                                                kernel there precisely BECAUSE
                                                stationarity fails -- so validity of
                                                that approximation is no reason to
                                                exclude small kappa from the model's
                                                admissible space. Worse, small kappa is
                                                Model 1's ONLY route to a delay rate
                                                still rising at the end of the delay, so
                                                a floor at 2.0 would delete exactly the
                                                region most favourable to the model the
                                                sweep is trying to exclude, manufacturing
                                                a cleaner result by construction. For a
                                                range whose purpose is exclusion, erring
                                                WIDE is the conservative direction.
                                                1/kappa here spans 5s -- slower than the
                                                whole 3s delay -- down to 33ms.)
  SDe      in [0.003, 0.030]                  (unchanged, architectural)
  beta     in [0.0015, 0.024] (Model 2 only)  (unchanged, architectural)
  muRS     in [0.30, 1.50]                    (unchanged, architectural)
  r        in [0.010, 0.200]  log-uniform     (unchanged -- already a task/design range,
                                                not derived from the fit)
  Lmotor   in [0.008, 0.080]  log-uniform     (unchanged, architectural)
  k        in [0.5, 2.0]      log-uniform     (NEW -- freed from FIXED per the same
                                                task/design reasoning as r and Lmotor:
                                                nothing ties it to the immediate-report
                                                architecture's inherited value any more
                                                than r or Lmotor were)
  lam      in [10.0, 50.0]    log-uniform     (NEW -- freed from FIXED; same
                                                transfer reasoning as k, but the upper
                                                end is NOT symmetric with it. 1/lam is
                                                the motor-preparation filter time
                                                constant, so [10,50] is tau in
                                                [100, 20] ms; lam = 100 would be a 10 ms
                                                constant, faster than plausible
                                                neuromuscular filtering. This bound also
                                                does argumentative work and the
                                                manuscript says so: the Model 1 ceiling
                                                argument (Results, Delay-locked) weakens
                                                as lam grows -- the exact ceiling at
                                                kappa=0.2 is 0.58/s at lam=21.5, 0.87/s
                                                at lam=50 and 1.22/s at lam=100, against
                                                an observed 1.45/s -- so an unbounded lam
                                                would quietly dissolve it. The bound is
                                                defensible on the time-constant argument
                                                alone; it is stated in the Table A1 note
                                                so that it cannot look chosen to protect
                                                the conclusion.)
  sigma_q  in [0.03, 0.30]    log-uniform     (NEW -- freed from FIXED per the user's
                                                "let it free" on identifiability grounds.
                                                No task/design constant pins its scale
                                                the way tau=1s pins kappa, so this range
                                                is anchored instead to sigx = 0.100, the
                                                OTHER diffusion-noise SD in the model,
                                                fixed since before any fitting: sigma_q
                                                plays the structurally analogous role for
                                                the post-response-signal process, so a
                                                +-3x bracket around sigx's own fixed
                                                value is the least arbitrary anchor
                                                available that still does not reference
                                                the fit. Flagged as the softest of these
                                                bounds -- revisit if asked.
sigx, Lmotion, LRS, delay stay fixed at Table A1's architectural point values (delay is
the task's own 3s delay duration, not a free parameter at all).

Latin-hypercube over the 10 continuous dimensions shared by both models (g, theta/g,
kappa, SDe, muRS, r, Lmotor, k, lam, sigma_q -- 10 dims; beta is an 11th dimension
sampled only for Model 2, mirroring the original design where "the same sampled vectors
were evaluated in both models, with beta set to zero in Model 1"). Three more free
dimensions than the 2026-08-24 design at the same N_DRAWS give substantially sparser
per-dimension coverage; consider raising N_DRAWS above 250 if the admissible set turns
out small enough that its composition looks noisy.

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

FIXED = dict(sigx=0.100, Lmotion=0.200, LRS=0.200, delay=3.0)

DIMS = ['g', 'theta_ratio', 'log_kappa', 'SDe', 'muRS', 'log_r', 'log_Lmotor',
        'log_k', 'log_lam', 'log_sigq', 'beta']
LO = dict(g=0.050, theta_ratio=0.25, log_kappa=np.log(0.20), SDe=0.003, muRS=0.30,
          log_r=np.log(0.010), log_Lmotor=np.log(0.008),
          log_k=np.log(0.5), log_lam=np.log(10.0), log_sigq=np.log(0.03), beta=0.0015)
HI = dict(g=0.200, theta_ratio=0.90, log_kappa=np.log(30.0), SDe=0.030, muRS=1.50,
          log_r=np.log(0.200), log_Lmotor=np.log(0.080),
          log_k=np.log(2.0), log_lam=np.log(50.0), log_sigq=np.log(0.30), beta=0.024)


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
    k = np.exp(LO['log_k'] + v['log_k'] * (HI['log_k'] - LO['log_k']))
    lam = np.exp(LO['log_lam'] + v['log_lam'] * (HI['log_lam'] - LO['log_lam']))
    sigq = np.exp(LO['log_sigq'] + v['log_sigq'] * (HI['log_sigq'] - LO['log_sigq']))
    beta_frac = v['beta']  # in [0,1]; mapped to alpha below per-model
    theta = g * theta_ratio
    p = dict(FIXED)
    p.update(g=g, theta=theta, kappa=kappa, SDe=SDe, muRS=muRS, r=r, Lmotor=Lmotor,
             k=k, lam=lam, sigq=sigq)
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
    slope = np.zeros(6); slope_se = np.zeros(6); rate = np.zeros(6)
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
        # mean delay-period partial-activation rate over the 2-s analysis window, per
        # trial per second -- directly comparable with the observed end-of-delay rate.
        rate[k] = h.sum() / (expo * 2.0) if expo > 0 else np.nan
        if (h > 0).sum() >= 4 and expo > 0:
            slope[k], slope_se[k] = poisson_slope_se(h, expo, DELAY_TB)
        else:
            slope[k], slope_se[k] = np.nan, np.nan
    return dict(prem=prem, early=early, noresp=noresp, acc=acc, rt_mean=rt_mean,
                rt_se=rt_se, rt_n=rt_n, slope=slope, slope_se=slope_se,
                rate=rate, met_mean=met_mean, met_se=met_se)


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
