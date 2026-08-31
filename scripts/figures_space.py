"""figures_space.py -- the sampled parameter space, with no fitted parameter anywhere.

  fig_ramp_vs_rate.pdf   the central question in one panel: across the admissible space,
                         can Model 1 produce a delay-period ramp AND activity at the
                         observed rate?  Every admissible draw is a point; the best
                         Model 1 candidates re-simulated at NTR=40000 are overlaid with
                         their Monte-Carlo error bars.
  fig_traj_median.pdf    the manuscript's Figure 1 recomputed at the MEDIAN of the
                         sampled ranges rather than at fitted parameters, with two
                         distinct bands: across trials at that point, and across the
                         admissible parameter space.
"""
import numpy as np, json, os, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS
import fused_sweep as F

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
C_M1, C_M2 = '#c1272d', '#3182bd'
OBS_RATE, OBS_LO, OBS_HI = 1.446, 0.659, 1.029
recs = [json.loads(l) for l in open('results/sweep_a.jsonl')]

# ---------------------------------------------------------------- 1. ramp vs rate
fig, ax = plt.subplots(figsize=(6.0, 4.4))
ax.axhspan(OBS_LO, OBS_HI, color='0.85', zorder=0)
ax.axvline(OBS_RATE, color='0.6', lw=.8, ls=':', zorder=1)
for tag, key, col, mk in [('Model 2', 'model2', C_M2, '^'), ('Model 1', 'model1', C_M1, 's')]:
    xs, ys = [], []
    for r in recs:
        m = r[key]
        if not m['admissible']:
            continue
        sl = np.asarray(m['raw']['slope'], float); rt = np.asarray(m['raw']['rate'], float)
        if np.isfinite(sl).any() and np.isfinite(rt).any():
            xs.append(np.nanmean(rt)); ys.append(np.nanmean(sl))
    ax.plot(xs, ys, mk, color=col, ms=3.2, alpha=.45, lw=0, label=f'{tag}, admissible draws')

RC = [json.loads(l) for l in open('results/ramp_check.jsonl')]
by = {}
for r in RC:
    by.setdefault(r['draw'], []).append(r)
for d, v in by.items():
    sl = np.array([x['slope'] for x in v]); rt = np.mean([x['rate'] for x in v])
    se = sl.std(ddof=1)/np.sqrt(sl.size) if sl.size > 1 else 0.0
    ax.errorbar(rt, sl.mean(), yerr=1.96*se, fmt='o', color='k', ms=5, capsize=2.5,
                zorder=5, mfc='w', mew=1.2,
                label='Model 1 best candidates,\nNTR=40000, 95% MC interval'
                      if d == list(by)[0] else None)
ax.set_xscale('symlog', linthresh=1e-3)
ax.set_xlabel('delay-period partial-activation rate (per trial per s)')
ax.set_ylabel('delay-period log-rate slope (s$^{-1}$)')
ax.set_title('Can Model 1 produce a ramp AND activity at the observed rate?', fontsize=9)
ax.text(1.55, OBS_HI, ' observed', fontsize=7, color='0.4', va='bottom')
ax.legend(frameon=False, loc='upper left', fontsize=6.8)
fig.savefig('plots/fig_ramp_vs_rate.pdf', bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- 2. median-point Fig 1
row = np.full(len(F.DIMS), 0.5)
Pmed, bfrac = F.to_param(row)
Delta = (F.FIXED['delay'] + F.FIXED['LRS']) - F.FIXED['Lmotion']
beta_med = np.exp(F.LO['log_beta'] + bfrac*(F.HI['log_beta'] - F.LO['log_beta']))
COLS = {0.0: '#4d4d4d', 0.12: '#7fa8c9', 0.50: '#1a5c8a'}
STY = {0.3: '--', 2.0: '-'}


def traj(P, model, alpha, n=600, seed=41):
    """mean gated signal per cell, on a common stimulus-offset-locked axis"""
    out = []
    for k, (c, D) in enumerate(CELLS):
        q = dict(P); q['alpha'] = alpha
        s = simulate(q, c, D, n, model=model, seed=seed + k, post=2.0)
        m = measures(s); ri = s['resp_i']; v = np.flatnonzero(m['valid'])
        Y = s['Z'] + P['g']
        t = np.arange(Y.shape[1])*DT - D                 # from stimulus offset
        M = np.array([np.where(np.arange(Y.shape[1]) <= ri[j], Y[j], np.nan) for j in v])
        out.append((t, np.nanmean(M, 0), np.nanstd(M, 0), c, D))
    return out


fig, axs = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True)
for r, (model, alpha, name, col) in enumerate(
        [(2, beta_med/Delta, 'Model 2 (with temporal preparation)', C_M2),
         (1, 0.0, 'Model 1 (maintenance only)', C_M1)]):
    ax = axs[r]
    # band across the admissible parameter space: 5-95% of the mean signal over draws
    key = 'model2' if model == 2 else 'model1'
    env = []
    for rr in recs[:60]:
        if not rr[key]['admissible']:
            continue
        tt, mu, sd, _, _ = traj(rr['params'], model,
                                rr['alpha'] if model == 2 else 0.0, n=120, seed=900)[4]
        env.append(np.interp(np.linspace(-0.3, 3.2, 200), tt, mu))
    if env:
        E = np.vstack(env)
        xx = np.linspace(-0.3, 3.2, 200)
        ax.fill_between(xx, np.nanpercentile(E, 5, 0), np.nanpercentile(E, 95, 0),
                        color=col, alpha=.13, lw=0,
                        label='5-95% across the admissible space')
    for tt, mu, sd, c, D in traj(Pmed, model, alpha):
        w = (tt >= -0.3) & (tt <= 3.2)
        ax.plot(tt[w], mu[w], STY[D], color=COLS[c], lw=1.1)
        ax.fill_between(tt[w], (mu-sd)[w], (mu+sd)[w], color=COLS[c], alpha=.06, lw=0)
    ax.axhline(Pmed['theta'], color='0.5', ls=':', lw=.8)
    ax.axhline(Pmed['g'], color='0.3', ls=':', lw=.8)
    ax.text(3.15, Pmed['theta'], r' $\theta$', fontsize=7, va='center', color='0.4')
    ax.text(3.15, Pmed['g'], ' g', fontsize=7, va='center', color='0.3')
    ax.set_ylabel('gated motor-preparatory signal')
    ax.set_title(name, fontsize=8.5, loc='left')
    ax.legend(frameon=False, loc='upper left', fontsize=6.5)
axs[1].set_xlabel('time from stimulus offset (s)')
fig.suptitle('At the median of the sampled ranges — no fitted parameter', fontsize=9.5)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig('plots/fig_traj_median.pdf', bbox_inches='tight'); plt.close(fig)
print('written: fig_ramp_vs_rate.pdf, fig_traj_median.pdf')
