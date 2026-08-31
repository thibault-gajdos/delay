"""
model1_bound.py -- how much log-rate ramp can Model 1 produce from the variance
transient alone?

Setting alpha = 0 and keeping the exact finite-rho moments, the rate depends on kappa
and on the distance to the gate and on nothing else: SD_eps cancels, because it scales
the distance and the standard deviation together.  Two observed quantities -- the rate at
the end of the delay window and its rise across the window -- therefore give two
equations in two unknowns, and the question of whether Model 1 can satisfy both is
decidable without any stationary approximation.

This replaces the earlier argument from the ceiling nu <= omega/2pi, which is the
alpha = 0 special case of the crossing rate and cannot be used to constrain kappa from
data generated with a ramp.
"""
import numpy as np
from scipy.stats import norm

LAM = 21.471                 # the Dendauw et al. (2024) value, Table A1's point value
LAM_RANGE = (10.0, 50.0)     # Table A1's sampled range once lambda is freed (2026-08-31)
T0, T1 = 1.05, 3.05          # elapsed maintenance time at the ends of the window
TARGET = 1.446               # observed rate at the last delay bin, cluster-corrected
OBSERVED_SLOPE = 0.700


def moments(t, kap, lam=None):
    LAM = LAM_ if (LAM_ := lam) is not None else globals()['LAM']
    sig = np.sqrt(2*kap); c = LAM*sig/(LAM-kap)
    V = c**2*((1-np.exp(-2*kap*t))/(2*kap)
              - 2*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
              + (1-np.exp(-2*LAM*t))/(2*LAM))
    Vd = c**2*(kap*(1-np.exp(-2*kap*t))/2
               - 2*kap*LAM*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
               + LAM*(1-np.exp(-2*LAM*t))/2)
    K = c*(np.exp(-kap*t)-np.exp(-LAM*t))
    return V, Vd, K**2/2


def rate(t, kap, d, lam=None):
    """non-stationary crossing rate at distance d, with alpha = 0"""
    V, Vd, C = moments(t, kap, lam)
    sV = np.sqrt(V)
    m = C/V*d
    s = np.sqrt(max(Vd - C**2/V, 1e-30))
    return norm.pdf(d/sV)/sV*(m*norm.cdf(m/s) + s*norm.pdf(m/s))


def best_slope(frac, kappas=None, lam=None):
    """largest log-rate slope Model 1 can produce while reaching frac x TARGET"""
    lam = LAM if lam is None else lam
    kappas = kappas if kappas is not None else np.r_[np.linspace(0.20, 3, 60),
                                                     np.linspace(3, 30, 40)]
    best, bk = -np.inf, None
    for kap in kappas:
        if abs(kap - lam) < 0.2:
            continue
        sV = np.sqrt(moments(T1, kap, lam)[0])
        ds = np.linspace(1e-4, 6*sV, 600)
        n1 = np.array([rate(T1, kap, d, lam) for d in ds])
        n0 = np.array([rate(T0, kap, d, lam) for d in ds])
        ok = n1 >= frac*TARGET
        if ok.any():
            v = (np.log(n1[ok]/np.maximum(n0[ok], 1e-300))/(T1-T0)).max()
            if v > best:
                best, bk = v, kap
    return best, bk


def ceiling(kap, lam=None, t=T0, nd=600):
    """exact (non-stationary) upper bound on Model 1's delay-period rate at (kappa, lam),
    maximized over the distance from the maintained level to the gate.  This is the
    quantity Results (Delay-locked) quotes; it is NOT the stationary sqrt(lam*kappa)/2pi
    of Corollary C.1, which understates it wherever the fluctuation has not relaxed."""
    sV = np.sqrt(moments(T1, kap, lam)[0])
    return max(rate(t, kap, d, lam) for d in np.linspace(1e-4, 6*sV, nd))


def worst_ceiling(kap=0.2, lams=None, n=25):
    """the ceiling at its worst (largest) over Table A1's sampled lambda range -- the
    form the manuscript states, quantified over the box rather than asserted at one
    lambda.  Returns (ceiling, lambda attaining it)."""
    lams = np.linspace(*LAM_RANGE, n) if lams is None else lams
    vals = [(ceiling(kap, L), L) for L in lams]
    return max(vals)


def relaxation(kap, rho=1.0):
    """fraction of the stationary variance reached after `rho` seconds of maintenance;
    this is the other fact Appendix B uses, and it is what licenses the stationary form
    of the crossing rate everywhere in the analysis window"""
    return moments(rho, kap)[0] / moments(1e4, kap)[0]


if __name__ == '__main__':
    print("Both numbers are quantified over Table A1's sampled lambda range "
          f"{LAM_RANGE}, not\nasserted at a single lambda (2026-08-31).\n")

    c, L = worst_ceiling(0.2)
    print("1. ceiling on Model 1's delay-period rate at the kappa it would need "
          "(kappa ~ 0.2)")
    print(f"{'lambda':>10} {'ceiling (1/s)':>15}")
    for lam in (10.0, 21.471, 30.0, 40.0, 50.0):
        print(f"{lam:10.3f} {ceiling(0.2, lam):15.3f}")
    print(f"   worst case over the range: {c:.3f} /s at lambda = {L:.1f}, "
          f"against an observed {TARGET} /s")
    print(f"   -> the manuscript states the ceiling as under 0.9 /s\n")

    print("2. largest log-rate slope Model 1 can produce, by the rate it reaches")
    print(f"{'fraction of obs rate':>21} {'max slope (1/s)':>16} {'at lambda':>10}")
    for f in (0.5, 0.25, 0.10):
        vals = [(best_slope(f, lam=L)[0], L) for L in np.linspace(*LAM_RANGE, 5)]
        v, L = max(vals)
        print(f"{f:21.2f} {v:16.4f} {L:10.1f}")
    print(f"{'observed':>21} {OBSERVED_SLOPE:16.3f}")
    print("\nModel 1 reaches the observed slope only at about a tenth of the observed")
    print("rate; at a quarter of it the largest available slope is 0.067 /s, and at half")
    print("of it 0.0002 /s.  It buys the ramp only by giving up the rate.\n")

    print("fraction of the stationary variance reached one second after commitment")
    for kap in (0.2, 2.3, 12.8):
        print(f"  kappa = {kap:5.2f}: {100*relaxation(kap):8.4f}%")
    print("  (this is why the small-kappa bound above is computed from the exact kernel,")
    print("   not from the stationary form -- and why it is NOT a reason to exclude")
    print("   small kappa from the sampled space; see fused_sweep.py's docstring.)")
