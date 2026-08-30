"""
theory_checks.py -- numerical verification of the two appendices, against direct
integration of the stochastic differential equations at a 1024 Hz step (the EMG/EEG
native sampling rate; matches sim.py's DT).

Appendix A, the rate: the Kac-Rice expression is checked against a simulated count of
upcrossings, and the occupation identity Phi(-a) = nu * Dbar against the simulated time
above the level.  The Euler bias on crossing counts is upward and shrinks with the step;
it affects the simulated rate, not the direction of any effect.

Appendix A, the observable: the ratio of merged clusters to raw upcrossings under the
25 ms rule, which is the factor any comparison between the formula and a measured rate
must carry.

Appendix B is verified by model1_bound.py, which bounds the log-rate slope Model 1 can
produce as a function of the rate it reaches, without any stationary approximation.

The residual condition dependence of motor execution time, from the filter transient
after the response signal, which is what makes invariance a consequence rather than an
assumption.
"""
import numpy as np
from scipy.stats import norm

LAM, DT = 21.471, 1.0 / 1024.0


def check_rate(kap, sde=0.03, a=1.5, T=400.0, n=40, merge_ms=25):
    """Simulate the stationary fluctuation, count upcrossings and clusters."""
    sV = sde * np.sqrt(LAM / (LAM + kap))
    g = a * sV
    se = sde * np.sqrt(2 * kap)
    N = int(T / DT)
    rng = np.random.default_rng(1)
    raw = mrg = 0
    occ = 0.0
    gap = int(merge_ms / (DT * 1000))
    for _ in range(n):
        e = rng.normal(0, sde); v = rng.normal(0, sV)
        z = rng.standard_normal(N) * np.sqrt(DT) * se
        Y = np.empty(N)
        for i in range(N):
            e += -kap * e * DT + z[i]
            v += LAM * (e - v) * DT
            Y[i] = v
        act = Y > g
        occ += act.mean()
        d = np.diff(act.astype(np.int8))
        st = np.flatnonzero(d == 1) + 1
        en = np.flatnonzero(d == -1) + 1
        raw += st.size
        if st.size:
            k = 1
            for j in range(1, st.size):
                prev = en[en < st[j]]
                if prev.size and st[j] - prev[-1] >= gap:
                    k += 1
            mrg += k
    om = np.sqrt(LAM * kap)
    return dict(rice=om / (2 * np.pi) * np.exp(-a**2 / 2), raw=raw / (n * T),
                merged=mrg / (n * T), occ=occ / n, phi=norm.cdf(-a))


def check_met(lam=21.93, mu=0.605, r=0.065):
    """Residual dependence of motor execution time on the time of activation onset."""
    f = lambda u: mu * (u - (1 - np.exp(-lam * u)) / lam)
    out = []
    for uc in (0.122, 0.148, 0.170, 0.198):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(uc + mid) - f(uc) < r:
                lo = mid
            else:
                hi = mid
        out.append((uc, lo))
    return out


if __name__ == '__main__':
    print("Appendix A: crossing rate and occupation, a = 1.5")
    print(f"{'kappa':>6} {'Rice':>9} {'simulated':>10} {'bias':>7} "
          f"{'clusters/raw':>13} {'Phi(-a)':>9} {'occupancy':>10}")
    for kap in (2.0, 5.0, 9.4):
        c = check_rate(kap)
        print(f"{kap:6.2f} {c['rice']:9.4f} {c['raw']:10.4f} "
              f"{100*(c['raw']/c['rice']-1):6.1f}% {c['merged']/c['raw']:13.3f} "
              f"{c['phi']:9.4f} {c['occ']:10.4f}")
    print("\nmotor execution time, deterministic filter transient")
    print(f"{'PMET - L_RS':>12} {'MET':>9}")
    for uc, met in check_met():
        print(f"{uc:12.3f} {1000*met:8.1f} ms")
    v = [m for _, m in check_met()]
    print(f"  spread over the observed range: {1000*(max(v)-min(v)):.1f} ms")
