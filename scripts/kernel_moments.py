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



def moments(t, kap, lam):
    """(Var V, Var Vdot, Cov^2/2) at elapsed maintenance time t, from the kernel."""
    LAM = lam
    sig = np.sqrt(2*kap); c = LAM*sig/(LAM-kap)
    V = c**2*((1-np.exp(-2*kap*t))/(2*kap)
              - 2*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
              + (1-np.exp(-2*LAM*t))/(2*LAM))
    Vd = c**2*(kap*(1-np.exp(-2*kap*t))/2
               - 2*kap*LAM*(1-np.exp(-(LAM+kap)*t))/(LAM+kap)
               + LAM*(1-np.exp(-2*LAM*t))/2)
    K = c*(np.exp(-kap*t)-np.exp(-LAM*t))
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
