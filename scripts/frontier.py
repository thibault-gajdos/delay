"""
frontier.py -- the ramp-rate frontier of BOTH models, computed exactly from the kernel.

Corollary C.5 of the manuscript: with alpha = 0 the delay-period log-rate slope is
  d(ln nu)/dtau = 0.5*(u'/u - s'/s) + (a^2/2)*(s'/s),   s = Var V, u = Var Vdot,
whose first term is NEGATIVE throughout the analysis window (u is governed by the fast
filter lambda and has saturated; s carries the slow mode kappa and has not).  The ramp
therefore lives entirely in the second term, which needs a large a, while the rate carries
exp(-a^2/2), which needs a small one.  The frontier is the upper envelope of the
achievable (rate, slope) pairs over the ranges of Table A1.

Computing it as an envelope rather than as a constrained optimisation at each rate is both
faster and more useful: the whole curve comes from one scatter over (kappa, lambda, a),
and the same scatter gives the (kappa, lambda, a) that attains each point.

Model 2 is computed the same way with the barrier moving: temporal preparation carries the
gate-facing signal towards the gate at rate alpha, which enters the Kac-Rice formula only
through the relative derivative Vdot + alpha, so `rate(..., alpha=)` handles it. The extra
dimension is the ramp height over the analysis window, in units of sigma_V; it is capped at
3.0, which is NOT an arbitrary choice -- it is the 75th percentile of beta/sigma_V over the
Model 2 draws the admissibility screen actually accepts (median 1.0, maximum 31.2). The cap
matters because the analytic calculation knows nothing about admissibility: without it,
Model 2's frontier would grow without bound simply by making the ramp larger, whereas in
the model a ramp large enough carries the signal past the response criterion during the
delay and shows up as premature responses, which the screen rejects.

The asymmetry between the two frontiers is the point, and it is worth stating plainly. For
Model 1 the frontier is a genuine bound: alpha = 0 leaves only kappa, lambda and a, all
three are swept, and nothing outside the calculation can push it up. For Model 2 it is an
upper envelope under a cap read off the screened sample, so it says what the architecture
can reach, not what it typically does.

Writes results/frontier.json: both envelopes, plus the numbers the manuscript quotes
(the largest slope at the observed rate, and the largest rate at the observed slope).
"""
import numpy as np, json, sys, os, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel_moments import moments, rate

T0, T1 = 0.8, 2.8          # ends of the delay analysis window, in elapsed maintenance time
# The slope is estimated the way the DATA's slope is estimated: an offset-Poisson GLM over
# the same twenty 0.1 s bins, not a two-point secant between the ends of the window. The
# distinction is not cosmetic for Model 1, whose log-rate is strongly concave: the fit
# weights the late, event-rich bins, where the transient has already flattened, and returns
# 45-75% of the secant. Comparing a secant against a fitted slope would credit Model 1 with
# a ramp its own data analysis would not have measured.
TB = np.arange(20) * 0.1 + 0.05 + T0
XB = TB - TB.mean()


def glm_slope(mu):
    """Offset-Poisson GLM slope of expected counts mu over the twenty bins."""
    X = np.column_stack([np.ones_like(XB), XB])
    b = np.array([np.log(max(mu.mean(), 1e-300)), 0.0])
    for _ in range(200):
        m = np.exp(np.clip(X @ b, -30, 30))
        z = X @ b + (mu - m) / np.maximum(m, 1e-12)
        A = X.T @ (m[:, None] * X) + 1e-12 * np.eye(2)
        bn = np.linalg.solve(A, X.T @ (m * z))
        if np.max(np.abs(bn - b)) < 1e-12:
            return bn[1]
        b = bn
    return b[1]
# The observed slope to compare a single homogeneous simulated process against is the MEAN
# OF THE PER-CELL, PER-PARTICIPANT slopes (0.892 /s; the six cell means span 0.659-1.029),
# not the slope of the profile pooled across participants (0.817). Pooling mixes processes
# with different distances to the gate, and the log of a mixture of exponentials is flatter
# than the mixture of the logs, so the pooled figure understates the slope of any one
# participant's process -- using it would make the frontier comparison look kinder to
# Model 1 than it should.
OBS_RATE, OBS_SLOPE = 0.709, 0.892
KAPPAS = np.r_[np.geomspace(0.20, 3.0, 110), np.geomspace(3.0, 30.0, 30)]
LAMS = np.geomspace(10.0, 50.0, 9)


RAMP_CAP = 3.0        # ramp height over the window, in sigma_V; see the docstring


def scatter(model=1, lams=None):
    """Achievable (rate at T1, log-rate slope) pairs over the sampled ranges."""
    R, S, P = [], [], []
    ramps = np.array([0.0]) if model == 1 else np.linspace(0.0, RAMP_CAP, 40)
    a_grid = np.linspace(0.02, 4.5, 260) if model == 1 else np.linspace(0.02, 3.5, 60)
    for lam in (LAMS if lams is None else lams):
        for kap in KAPPAS:
            if abs(kap - lam) < 0.3:
                continue
            sV = np.sqrt(moments(T1, kap, lam)[0])
            for a1 in a_grid:
                d1 = a1 * sV                       # distance at the END of the window
                for b in ramps:
                    al = b * sV / (T1 - T0)
                    mu = np.array([rate(t, kap, d1 + al * (T1 - t), lam, al) for t in TB])
                    if not np.all(mu > 0):
                        continue
                    R.append(mu[-1]); S.append(glm_slope(mu))
                    P.append((kap, lam, a1, b))
    return np.array(R), np.array(S), P


def envelope(R, S, P, nbin=90):
    lo, hi = np.log(R.min()), np.log(R.max())
    edges = np.linspace(lo, hi, nbin + 1)
    idx = np.digitize(np.log(R), edges) - 1
    out = []
    for b in range(nbin):
        m = np.flatnonzero(idx == b)
        if m.size == 0:
            continue
        j = m[int(np.argmax(S[m]))]
        out.append(dict(rate=float(R[j]), slope=float(S[j]), kappa=float(P[j][0]),
                        lam=float(P[j][1]), a=float(P[j][2]), ramp=float(P[j][3])))
    out.sort(key=lambda r: r['rate'])
    # make it a true frontier: the largest slope achievable at a rate of AT LEAST x is
    # non-increasing in x, so take the running maximum from the high-rate end
    best = -np.inf
    for r in reversed(out):
        best = max(best, r['slope']); r['slope_at_least'] = float(best)
    return out


def _queries(R, S, P):
    def at_rate(f):
        m = R >= f
        if not m.any():
            return None
        j = int(np.argmax(np.where(m, S, -np.inf)))
        return dict(rate=float(R[j]), slope=float(S[j]), kappa=float(P[j][0]),
                    lam=float(P[j][1]), a=float(P[j][2]), ramp=float(P[j][3]))

    def at_slope(g):
        m = S >= g
        if not m.any():
            return None
        j = int(np.argmax(np.where(m, R, -np.inf)))
        return dict(rate=float(R[j]), slope=float(S[j]), kappa=float(P[j][0]),
                    lam=float(P[j][1]), a=float(P[j][2]), ramp=float(P[j][3]))
    return at_rate, at_slope


FLOORS = [0.085, 0.2, 0.3, 0.4, 0.5, 0.709, 1.0]

if __name__ == '__main__':
    # Sharded over the lambda grid so the scatter can be split across cores:
    #   FR_MODEL=1 FR_SHARD_N=9 FR_SHARD_I=$i python3 frontier.py  -> results/frontier_m1_sh$i.npz
    # then `python3 frontier.py merge` combines them into results/frontier.json. Sharding by
    # lambda is exact: the scatter is a product grid, no shard needs another's points, and the
    # envelope of a union is the envelope of the per-shard point sets concatenated.
    if len(sys.argv) > 1 and sys.argv[1] == 'merge':
        res = dict(window=[T0, T1], obs_rate=OBS_RATE, obs_slope=OBS_SLOPE,
                   ramp_cap_sigmaV=RAMP_CAP,
                   slope_statistic='offset-Poisson GLM over the same twenty 0.1 s bins the '
                                   'data analysis uses, not a two-point secant',
                   models={})
        for model in (1, 2):
            fs = sorted(glob.glob(f'results/frontier_m{model}_sh*.npz'))
            if not fs:
                print(f'--- Model {model}: no shards, skipped')
                continue
            R = np.concatenate([np.load(f)['R'] for f in fs])
            S = np.concatenate([np.load(f)['S'] for f in fs])
            P = [tuple(x) for x in np.concatenate([np.load(f)['P'] for f in fs])]
            at_rate, at_slope = _queries(R, S, P)
            res['models'][f'model{model}'] = dict(
                n_points=int(R.size), n_shards=len(fs), envelope=envelope(R, S, P),
                max_slope_at_observed_rate=at_rate(OBS_RATE),
                max_rate_at_observed_slope=at_slope(OBS_SLOPE),
                max_slope_by_rate_floor={str(f): at_rate(f) for f in FLOORS})
            print(f'--- Model {model}: {R.size} points from {len(fs)} shards')
            for f in FLOORS:
                v = at_rate(f)
                if v:
                    print('  rate >= %5.3f : max slope %+8.3f  (kappa=%.2f, lam=%.0f, a=%.2f, ramp=%.2f)'
                          % (f, v['slope'], v['kappa'], v['lam'], v['a'], v['ramp']))
            v = at_slope(OBS_SLOPE)
            print('  largest rate at the observed slope %+.3f : %s'
                  % (OBS_SLOPE, ('%.3f/s' % v['rate']) if v else 'unattainable anywhere'))
        json.dump(res, open('results/frontier.json', 'w'), indent=1)
        print('\nwrote results/frontier.json')
        sys.exit(0)

    model = int(os.environ.get('FR_MODEL', '1'))
    n = int(os.environ.get('FR_SHARD_N', '1'))
    i = int(os.environ.get('FR_SHARD_I', '0'))
    lams = LAMS[i::n]
    t0 = time.time()
    R, S, P = scatter(model, lams=lams)
    out = f'results/frontier_m{model}_sh{i}.npz'
    np.savez(out, R=R, S=S, P=np.array(P, float))
    print(f'model {model} shard {i}/{n}: {R.size} points over lambda={np.round(lams, 1)} '
          f'in {time.time() - t0:.0f}s -> {out}')
