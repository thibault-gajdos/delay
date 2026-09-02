"""
kernel_moments.py -- the exact (non-stationary) second moments of the maintained
fluctuation V, and the crossing rate that follows from them.

V is the output of the first-order filter (rate lambda) driven by the Ornstein-Uhlenbeck
process (rate kappa), started from rest at the moment maintenance begins: V(tau0) = 0 with
zero variance. Solving the two linear stochastic differential equations in sequence gives V
as a Wiener integral against a kernel, from which the variance of V, the variance of its
derivative and their covariance follow in closed form. `moments` returns those three;
`rate` is the Kac-Rice upcrossing rate at a fixed distance d, keeping the V/Vdot covariance
that non-stationarity induces.

These are computed from the kernel rather than from the stationary limits BECAUSE the
process is not stationary over most of the sampled range of kappa: the relative departure
from the stationary values one second into maintenance is of order exp(-2 kappa), below
1e-9 only above kappa ~ 11, whereas Table A1 samples down to kappa = 0.2. Small kappa is
also the region most favourable to Model 1, so it is exactly the region a bound must not
approximate away.

Used by frontier.py (the ramp-rate frontier) and by theory_checks.py.
"""
import numpy as np
from scipy.stats import norm



DEGENERATE = 1e-5      # relative |lam - kap| at which to switch to the repeated-pole forms


def moments(t, kap, lam):
    """(Var H, Var Hdot, Cov(H, Hdot)) at elapsed maintenance time t, from the kernel.

    The kernel is the impulse response of the two filters in series, K(v) = lam sig
    (e^{-kap v} - e^{-lam v})/(lam - kap), and the three moments are integrals of it:
    Var H = int K^2, Var Hdot = int K'^2, Cov = int K K' = K(t)^2/2.

    kap = lam is a REPEATED POLE, not a singularity of the process: the two modes of the
    cascade coincide, the basis {e^{-kap v}, e^{-lam v}} degenerates into
    {e^{-lam v}, v e^{-lam v}}, and the kernel becomes K(v) = lam sig v e^{-lam v}. The
    limit exists; only the parameterisation by two distinct modes breaks down.

    The general expressions are 0/0 there and lose precision to cancellation as they
    approach it, while the repeated-pole forms lose accuracy as they move away from it.
    Measured against numerical quadrature of the kernel, the two error curves cross near a
    relative separation of 3e-6: at 1e-3 the general form is exact to 2e-10 and the
    degenerate one wrong by 1.5e-3; at 1e-7 the general form returns 5.55 instead of 0.5
    while the degenerate one is exact to 1.5e-7. Switching at 1e-5 keeps the worst-case
    relative error of either branch below 1.5e-5.
    """
    sig = np.sqrt(2*kap)
    if abs(lam - kap) <= DEGENERATE*lam:
        E = np.exp(-2*lam*t)
        V = sig**2/(4*lam)*(1 - E*(1 + 2*lam*t + 2*lam**2*t**2))
        Vd = lam*sig**2/4  *(1 - E*(1 - 2*lam*t + 2*lam**2*t**2))
        K = lam*sig*t*np.exp(-lam*t)
        return V, Vd, K**2/2
    c = lam*sig/(lam-kap)
    V = c**2*((1-np.exp(-2*kap*t))/(2*kap)
              - 2*(1-np.exp(-(lam+kap)*t))/(lam+kap)
              + (1-np.exp(-2*lam*t))/(2*lam))
    Vd = c**2*(kap*(1-np.exp(-2*kap*t))/2
               - 2*kap*lam*(1-np.exp(-(lam+kap)*t))/(lam+kap)
               + lam*(1-np.exp(-2*lam*t))/2)
    K = c*(np.exp(-kap*t)-np.exp(-lam*t))
    return V, Vd, K**2/2


def rate(t, kap, d, lam, alpha=0.0):
    """Upcrossing rate of the barrier at distance d, at elapsed maintenance time t.

    `alpha` is the rate at which temporal preparation moves the signal towards the gate,
    in the same units as d per second (Model 2; alpha = 0 is Model 1). A moving barrier
    enters the Kac-Rice formula through the relative derivative Vdot - ddot = Vdot + alpha,
    so it shifts the conditional mean of that derivative by +alpha and nothing else; the
    density at the barrier is unchanged. `d` must already be the distance AT TIME t, i.e.
    d0 - alpha*t, not the initial one.
    """
    V, Vd, C = moments(t, kap, lam)
    sV = np.sqrt(V)
    m = C/V*d + alpha
    s = np.sqrt(max(Vd - C**2/V, 1e-30))
    return norm.pdf(d/sV)/sV*(m*norm.cdf(m/s) + s*norm.pdf(m/s))
