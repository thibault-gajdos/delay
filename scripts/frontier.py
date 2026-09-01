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
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel_moments import moments, rate

T0, T1 = 0.8, 2.8          # ends of the delay analysis window, in elapsed maintenance time
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


def scatter(model=1):
    """Achievable (rate at T1, log-rate slope) pairs over the sampled ranges."""
    R, S, P = [], [], []
    ramps = np.array([0.0]) if model == 1 else np.linspace(0.0, RAMP_CAP, 40)
    a_grid = np.linspace(0.02, 4.5, 260) if model == 1 else np.linspace(0.02, 3.5, 60)
    for lam in LAMS:
        for kap in KAPPAS:
            if abs(kap - lam) < 0.3:
                continue
            sV = np.sqrt(moments(T1, kap, lam)[0])
            for a1 in a_grid:
                d1 = a1 * sV                       # distance at the END of the window
                for b in ramps:
                    al = b * sV / (T1 - T0)
                    n1 = rate(T1, kap, d1, lam, al)
                    n0 = rate(T0, kap, d1 + al * (T1 - T0), lam, al)
                    if n0 <= 0 or n1 <= 0:
                        continue
                    R.append(n1); S.append(np.log(n1 / n0) / (T1 - T0))
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


if __name__ == '__main__':
    res = dict(window=[T0, T1], obs_rate=OBS_RATE, obs_slope=OBS_SLOPE,
               ramp_cap_sigmaV=RAMP_CAP, models={})
    for model in (1, 2):
        R, S, P = scatter(model)
        env = envelope(R, S, P)

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

        floors = [0.085, 0.2, 0.3, 0.4, 0.5, 0.709, 1.0]
        res['models'][f'model{model}'] = dict(
            n_points=int(R.size), envelope=env,
            max_slope_at_observed_rate=at_rate(OBS_RATE),
            max_rate_at_observed_slope=at_slope(OBS_SLOPE),
            max_slope_by_rate_floor={str(f): at_rate(f) for f in floors})
        print(f"--- Model {model}: {R.size} points")
        for f in floors:
            v = at_rate(f)
            if v:
                print("  rate >= %5.3f : max slope %+8.3f  (kappa=%.2f, lam=%.0f, a=%.2f, ramp=%.2f)"
                      % (f, v['slope'], v['kappa'], v['lam'], v['a'], v['ramp']))
        v = at_slope(OBS_SLOPE)
        if v:
            print("  largest rate at the observed slope %+.3f : %.3f/s" % (OBS_SLOPE, v['rate']))
    json.dump(res, open('results/frontier.json', 'w'), indent=1)
    print('\nwrote results/frontier.json')
