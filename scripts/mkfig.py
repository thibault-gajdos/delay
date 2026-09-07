"""Appendix figure: how far the stationary rate is from the exact rate,
over the sampled parameter domain.  Analytic -- no simulation needed."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'axes.spines.top': False, 'axes.spines.right': False})
Phi, phi = norm.cdf, norm.pdf
C_M1, C_M2 = '#c1272d', '#3182bd'

def moments(r, l, k):
    """Warm-start second-order moments of H:  sigma_H^2, c_H, sigma_Hdot^2."""
    s2 = 1/(2*k); Ai = l*s2/(k+l); Bi = l**2*s2*k/(k+l)
    E1, E2 = np.exp(-(l+k)*r), np.exp(-2*l*r)
    X = (l+k)*(E2-E1)/(l-k) if abs(l-k) > 1e-12 else -2*l*r*E2
    return Ai*(1+X-E1), -l*X*Ai, Bi*(1+(l/k)*(X+E1)), Ai, Bi

def rate(r, l, k, bd, d):
    A, cH, B, _, _ = moments(r, l, k); sH = np.sqrt(A)
    M = bd + cH/A*d; S = np.where(B - cH**2/A > 0, B - cH**2/A, np.nan)
    u = M/np.sqrt(S)
    return phi(d/sH)/sH*(M*Phi(u) + np.sqrt(S)*phi(u))

def rate_stat(l, k, bd, d):
    s2 = 1/(2*k); Ai = l*s2/(k+l); Bi = l**2*s2*k/(k+l); us = bd/np.sqrt(Bi)
    return phi(d/np.sqrt(Ai))/np.sqrt(Ai)*(bd*Phi(us) + np.sqrt(Bi)*phi(us))

rho = np.geomspace(0.02, 3.0, 500)
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
for ax, mod, col, nm in [(axes[0], 'model1', C_M1, 'Model 1'),
                         (axes[1], 'model2', C_M2, 'Model 2')]:
    n, worst = 0, 0.0
    for line in open('results/sweep_screened.jsonl'):
        j = json.loads(line); p = j['params']
        l, k = p['lam'], p['kappa']
        beta = j['beta'] if mod == 'model2' else 0.0
        bd = beta/(3.0 + p['LRS'] - p['Lmotion'])
        d = (p['g'] - p['theta']) - bd*rho
        if np.any(d <= 0): continue
        with np.errstate(all='ignore'):
            e = np.abs(rate(rho, l, k, bd, d)/rate_stat(l, k, bd, d) - 1)
        if not np.all(np.isfinite(e)): continue
        ax.plot(rho, np.maximum(100*e, 1e-7), color=col, lw=.35, alpha=.15, zorder=1)
        n += 1; worst = max(worst, e[rho >= 1.0].max())
    ax.axvspan(1.0, 3.0, color='0.90', zorder=0)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(0.02, 3); ax.set_ylim(1e-6, 300)
    ax.axhline(1.0, color='0.55', lw=.7, ls=':', zorder=2)
    ax.set_xlabel('elapsed maintenance time (s)')
    ax.set_title(f'{nm}, {n} parameter draws', fontsize=8)
    ax.set_xticks([0.02, 0.1, 0.3, 1, 3]); ax.set_xticklabels(['0.02','0.1','0.3','1','3'])
    print(f'{nm}: n={n}, largest error inside the window = {100*worst:.2e} %')
axes[0].set_ylabel(r'$|\,\nu\,/\,\nu_{\rm stat}-1\,|$   (%)')
axes[0].text(1.06, 90, 'analysis\nwindow', fontsize=6.5, color='0.35', va='top')
axes[0].text(0.022, 1.6, '1%', fontsize=6.5, color='0.45')
fig.savefig('/tmp/approx_error.pdf', bbox_inches='tight')
fig.savefig('/tmp/approx_error.png', dpi=200, bbox_inches='tight')
