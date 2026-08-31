"""
frontier.py -- Model 1's ramp-rate frontier, computed exactly from the kernel.

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

Writes results/frontier.json: the envelope, plus the two numbers the manuscript quotes
(the largest slope at the observed rate, and the largest rate at the observed slope).
"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model1_bound import moments, rate

T0, T1 = 0.8, 2.8          # ends of the delay analysis window, in elapsed maintenance time
OBS_RATE, OBS_SLOPE = 0.709, 0.817
KAPPAS = np.r_[np.geomspace(0.20, 3.0, 110), np.geomspace(3.0, 30.0, 30)]
LAMS = np.geomspace(10.0, 50.0, 9)


def scatter():
    R, S, P = [], [], []
    for lam in LAMS:
        for kap in KAPPAS:
            if abs(kap - lam) < 0.3:
                continue
            sV = np.sqrt(moments(T1, kap, lam)[0])
            for a in np.linspace(0.02, 4.5, 260):
                d = a * sV
                n1 = rate(T1, kap, d, lam); n0 = rate(T0, kap, d, lam)
                if n0 <= 0 or n1 <= 0:
                    continue
                R.append(n1); S.append(np.log(n1 / n0) / (T1 - T0)); P.append((kap, lam, a))
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
        out.append(dict(rate=float(R[j]), slope=float(S[j]),
                        kappa=float(P[j][0]), lam=float(P[j][1]), a=float(P[j][2])))
    out.sort(key=lambda r: r['rate'])
    # make it a true frontier: the largest slope achievable at a rate of AT LEAST x is
    # non-increasing in x, so take the running maximum from the high-rate end
    best = -np.inf
    for r in reversed(out):
        best = max(best, r['slope']); r['slope_at_least'] = float(best)
    return out


if __name__ == '__main__':
    R, S, P = scatter()
    env = envelope(R, S, P)
    # queries answered from the raw scatter, not from the binned envelope: binning can
    # only lose the maximum, and this frontier is used to EXCLUDE a model, so it must not
    # be under-reported.
    def max_slope_at_least(f):
        m = R >= f
        if not m.any():
            return None
        j = int(np.argmax(np.where(m, S, -np.inf)))
        return dict(rate=float(R[j]), slope=float(S[j]), kappa=float(P[j][0]),
                    lam=float(P[j][1]), a=float(P[j][2]))

    def max_rate_at_slope(g):
        m = S >= g
        if not m.any():
            return None
        j = int(np.argmax(np.where(m, R, -np.inf)))
        return dict(rate=float(R[j]), slope=float(S[j]), kappa=float(P[j][0]),
                    lam=float(P[j][1]), a=float(P[j][2]))

    at_rate = max_slope_at_least(OBS_RATE)
    at_slope = max_rate_at_slope(OBS_SLOPE)
    table = {f'{f}': max_slope_at_least(f) for f in [0.085, 0.2, 0.3, 0.4, 0.5, 0.709, 1.0]}
    res = dict(window=[T0, T1], obs_rate=OBS_RATE, obs_slope=OBS_SLOPE,
               n_points=int(R.size), envelope=env,
               max_slope_at_observed_rate=at_rate,
               max_rate_at_observed_slope=at_slope,
               max_slope_by_rate_floor=table)
    json.dump(res, open('results/frontier.json', 'w'), indent=1)
    print(f"{R.size} (kappa,lambda,a) points; envelope has {len(env)} bins")
    print("largest slope at a rate of at least %.3f/s : %+.4f  (kappa=%.2f, lam=%.0f, a=%.2f)"
          % (OBS_RATE, at_rate['slope'], at_rate['kappa'], at_rate['lam'], at_rate['a']))
    if at_slope:
        print("largest rate at a slope of at least %+.3f  : %.4f/s (kappa=%.2f, lam=%.0f, a=%.2f)"
              % (OBS_SLOPE, at_slope['rate'], at_slope['kappa'], at_slope['lam'], at_slope['a']))
    for f, v in table.items():
        if v:
            print("  rate >= %5.3f : max slope %+.3f  (kappa=%.2f, lam=%.0f, a=%.2f)"
                  % (float(f), v['slope'], v['kappa'], v['lam'], v['a']))
