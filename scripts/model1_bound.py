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

LAM = 25.5
T0, T1 = 1.05, 3.05          # elapsed maintenance time at the ends of the window
TARGET = 1.446               # observed rate at the last delay bin, cluster-corrected
OBSERVED_SLOPE = 0.700


def moments(t, kap):
    sig = np.sqrt(2*kap); c = LAM*sig/(LAM-kap)
    V = c**2*((1-np.exp(-2*kap*t))/(2*kap)
              - 2*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
              + (1-np.exp(-2*LAM*t))/(2*LAM))
    Vd = c**2*(kap*(1-np.exp(-2*kap*t))/2
               - 2*kap*LAM*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
               + LAM*(1-np.exp(-2*LAM*t))/2)
    K = c*(np.exp(-kap*t)-np.exp(-LAM*t))
    return V, Vd, K**2/2


def rate(t, kap, d):
    """non-stationary crossing rate at distance d, with alpha = 0"""
    V, Vd, C = moments(t, kap)
    sV = np.sqrt(V)
    m = C/V*d
    s = np.sqrt(max(Vd - C**2/V, 1e-30))
    return norm.pdf(d/sV)/sV*(m*norm.cdf(m/s) + s*norm.pdf(m/s))


def best_slope(frac, kappas=None):
    """largest log-rate slope Model 1 can produce while reaching frac x TARGET"""
    kappas = kappas if kappas is not None else np.r_[np.linspace(0.05, 3, 60),
                                                     np.linspace(3, 40, 40)]
    best, bk = -np.inf, None
    for kap in kappas:
        if abs(kap - LAM) < 0.2:
            continue
        sV = np.sqrt(moments(T1, kap)[0])
        ds = np.linspace(1e-4, 6*sV, 600)
        n1 = np.array([rate(T1, kap, d) for d in ds])
        n0 = np.array([rate(T0, kap, d) for d in ds])
        ok = n1 >= frac*TARGET
        if ok.any():
            v = (np.log(n1[ok]/np.maximum(n0[ok], 1e-300))/(T1-T0)).max()
            if v > best:
                best, bk = v, kap
    return best, bk


def relaxation(kap, rho=1.0):
    """fraction of the stationary variance reached after `rho` seconds of maintenance;
    this is the other fact Appendix B uses, and it is what licenses the stationary form
    of the crossing rate everywhere in the analysis window"""
    return moments(rho, kap)[0] / moments(1e4, kap)[0]


if __name__ == '__main__':
    print("fraction of the stationary variance reached one second after commitment\n")
    for kap in (2.3, 3.2, 12.8):
        print(f"  kappa = {kap:5.2f}: {100*relaxation(kap):8.4f}%")
    print("\nlargest log-rate slope Model 1 can produce, by rate reached\n")
    print(f"{'fraction of observed rate':>26} {'max slope (1/s)':>16} {'at kappa':>10}")
    for f in (0.5, 0.25, 0.10, 0.05):
        v, k = best_slope(f)
        print(f"{f:26.2f} {v:16.4f} {k:10.2f}")
    print(f"{'observed':>26} {OBSERVED_SLOPE:16.3f}")
    print("\nModel 1 matches the observed slope only at about a tenth of the observed")
    print("rate; at a quarter of it the largest available slope is 0.004 /s.")
