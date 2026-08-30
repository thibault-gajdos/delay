"""figures.py -- figures for the note.  Run after calibrate.py has produced the fits."""
import numpy as np, json, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, rate_profile, DT, _bursts
from delaymodel.targets import CELLS, LAB, build, poisson_slope

COL = {0: '#9ecae1', 12: '#3182bd', 50: '#08306b'}
CC = [0, 12, 50]
STY = {0.3: '-', 2.0: '--'}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'lines.linewidth': 1.1})


def cell_col(k):
    return COL[CC[k % 3]]


def cell_sty(k):
    return '-' if k < 3 else '--'


# ---------------------------------------------------------------- trajectories
def fig_trajectories(P, model=2, n=800, fname='plots/fig_traj.pdf'):
    """Simulated counterpart of the manuscript's Figure 1: mean gated signal,
    time-locked to the four trial events, with a one standard deviation band."""
    fig, ax = plt.subplots(1, 4, figsize=(9.2, 2.7),
                           gridspec_kw=dict(width_ratios=[2.2, 3, 0.9, 0.9], wspace=0.12))
    for k, (c, D) in enumerate(CELLS):
        s = simulate(P, c, D, n, model, seed=41 + k)
        m = measures(s)
        Y = s['Z'] + P['g']                     # back to the gated signal
        v = m['valid']; ri = s['resp_i']
        col, sty = cell_col(k), cell_sty(k)
        t = np.arange(Y.shape[1]) * DT
        # mask each trial after its response
        M = np.array([np.where(np.arange(Y.shape[1]) <= (ri[j] if ri[j] >= 0 else 10**9),
                               Y[j], np.nan) for j in np.flatnonzero(v)])
        mu = np.nanmean(M, 0); sd = np.nanstd(M, 0)
        w = t <= 2.2
        ax[0].plot(t[w], mu[w], sty, color=col)
        ax[0].fill_between(t[w], mu[w] - sd[w], mu[w] + sd[w], color=col, alpha=.07, lw=0)
        t2 = t - D
        w2 = (t2 >= -0.1) & (t2 <= 3 + P['LRS'])
        ax[1].plot(t2[w2], mu[w2], sty, color=col)
        ax[1].fill_between(t2[w2], mu[w2] - sd[w2], mu[w2] + sd[w2], color=col, alpha=.07, lw=0)
        t3 = t - (D + 3)
        w3 = (t3 >= 0) & (t3 <= 0.45)
        ax[2].plot(t3[w3], mu[w3], sty, color=col)
        # response-locked
        L = int(0.25 / DT)
        R = np.full((v.sum(), L), np.nan)
        for i, j in enumerate(np.flatnonzero(v)):
            if ri[j] >= L:
                R[i] = Y[j, ri[j] - L:ri[j]]
        ax[3].plot(np.arange(-L, 0) * DT, np.nanmean(R, 0), sty, color=col)
    for a, ttl in zip(ax, ['from stimulus onset', 'from stimulus offset',
                           'from response signal', 'from response']):
        a.axhline(P['g'], color='0.3', ls=':', lw=.8)
        a.axhline(P['theta'], color='0.6', ls=':', lw=.8)
        a.set_title(ttl, fontsize=8)
        a.set_ylim(-0.01, P['g'] + P['r'] * 1.6)
    ax[0].set_ylabel('gated motor-preparatory signal')
    ax[0].text(0.02, P['g'] + 0.004, '$g$', fontsize=8, color='0.3')
    ax[0].text(0.02, P['theta'] + 0.004, r'$\theta$', fontsize=8, color='0.6')
    for a in ax[1:]:
        a.set_yticklabels([])
    h = [plt.Line2D([], [], color=COL[c], ls='-') for c in CC] + \
        [plt.Line2D([], [], color='0.4', ls=s) for s in ['-', '--']]
    ax[1].legend(h, ['0%', '12%', '50%', '0.3 s', '2 s'], ncol=5, fontsize=6.5,
                 frameon=False, loc='upper left')
    ax[0].set_xlabel('from onset (s)'); ax[1].set_xlabel('from offset (s)')
    ax[2].set_xlabel('(s)'); ax[3].set_xlabel('(s)')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- single trial
def fig_single(P, n=400, fname='plots/fig_single.pdf'):
    s = simulate(P, 0.12, 2.0, n, 2, seed=7)
    m = measures(s)
    cand = [j for j in np.flatnonzero(m['valid'])
            if len(m['events'][j]) >= 1 and any(e > s['D'] + 1 for e in m['events'][j])]
    j = cand[0]
    Y = s['Z'][j] + P['g']; t = np.arange(len(Y)) * DT
    ri = s['resp_i'][j]
    fig, ax = plt.subplots(2, 1, figsize=(7.4, 3.2), sharex=True,
                           gridspec_kw=dict(height_ratios=[3, 1], hspace=.12))
    ax[0].plot(t[:ri + 1], Y[:ri + 1], color='#08306b', lw=.7)
    ax[0].axhline(P['g'], color='0.3', ls=':', lw=.8)
    ax[0].axhline(P['g'] + P['r'], color='#c1272d', ls=':', lw=.8)
    ax[0].text(0.05, P['g'] + 0.002, '$g$', fontsize=8, color='0.3')
    ax[0].text(0.05, P['g'] + P['r'] + 0.002, '$g+r$', fontsize=8, color='#c1272d')
    for x, lab in [(s['D'], 'offset'), (s['D'] + 3, 'signal')]:
        for a in ax:
            a.axvline(x, color='0.7', lw=.7)
        ax[0].text(x + .02, ax[0].get_ylim()[1] * .95, lab, fontsize=6.5, color='0.5')
    act = Y[:ri + 1] > P['g']
    for a, b in _bursts(act):
        col = '#c1272d' if b >= ri - 1 else '#e6892f'
        ax[1].axvspan(a * DT, b * DT, color=col, alpha=.85)
    ax[1].set_yticks([]); ax[1].set_ylabel('detected', fontsize=7)
    ax[1].set_xlabel('time from stimulus onset (s)')
    ax[0].set_ylabel('gated signal')
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- rates
def fig_rates(P, E, n=900, fname='plots/fig_rates.pdf'):
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 2.7), sharey=True,
                           gridspec_kw=dict(width_ratios=[1, 1], wspace=.08))
    for k, (c, D) in enumerate(CELLS):
        s = simulate(P, c, D, n, 2, seed=131 + k); m = measures(s)
        md = rate_profile(m, D + 1, D + 3, 20)
        ms = rate_profile(m, 0.0, 2.0, 20)
        tb = np.arange(20) * 0.1 + 0.05
        col, sty = cell_col(k), cell_sty(k)
        ax[0].plot(tb, E['stim'].mean(0)[k], sty, color=col, alpha=.35, lw=.8)
        ax[0].plot(tb, ms, sty, color=col, lw=1.4)
        ax[1].plot(tb + 1, E['dly'].mean(0)[k], sty, color=col, alpha=.35, lw=.8)
        ax[1].plot(tb + 1, md, sty, color=col, lw=1.4)
    ax[0].set_title('stimulus period', fontsize=8)
    ax[1].set_title('delay period', fontsize=8)
    ax[0].set_ylabel('partial activations per 0.1 s bin')
    ax[0].set_xlabel('time from stimulus onset (s)')
    ax[1].set_xlabel('time from stimulus offset (s)')
    h = [plt.Line2D([], [], color='0.4', alpha=.35, lw=.8), plt.Line2D([], [], color='0.4', lw=1.4)]
    ax[1].legend(h, ['data', 'model'], fontsize=6.5, frameon=False)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- shape tests
def fig_shape(E, fname='plots/fig_shape.pdf'):
    tb = np.arange(20) * 0.1 + 0.05; tc = tb - tb.mean()
    lr = np.log(E['dly'].mean(0) .mean(-1) / 0.1)
    sl = np.nanmean(E['slope'], 0)
    fig, ax = plt.subplots(1, 3, figsize=(9.0, 2.5), gridspec_kw=dict(wspace=.30))
    for k in range(6):
        y = E['dly'].mean(0)[k]
        ax[0].plot(tb + 1, np.log(y / 0.1), cell_sty(k), color=cell_col(k))
    ax[0].set_xlabel('time from stimulus offset (s)'); ax[0].set_ylabel(r'$\ln\nu$')
    ax[0].set_title('delay log-rate profiles', fontsize=8)
    for k in range(6):
        ax[1].plot(lr[k], sl[k], 'o', color=cell_col(k),
                   mfc=cell_col(k) if k >= 3 else 'w', ms=5)
    b = np.polyfit(lr, sl, 1)
    xs = np.linspace(lr.min(), lr.max(), 10)
    ax[1].plot(xs, np.polyval(b, xs), '-', color='0.5', lw=.8)
    ax[1].set_xlabel(r'$\ln\nu$'); ax[1].set_ylabel('log-rate slope (s$^{-1}$)')
    ax[1].set_title(f'P2: slope falls as rate rises', fontsize=8)
    for k in range(6):
        ax[2].plot(sl[k] ** 2, lr[k], 'o', color=cell_col(k),
                   mfc=cell_col(k) if k >= 3 else 'w', ms=5)
    b2 = np.polyfit(sl ** 2, lr, 1)
    xs = np.linspace((sl ** 2).min(), (sl ** 2).max(), 10)
    ax[2].plot(xs, np.polyval(b2, xs), '-', color='0.5', lw=.8)
    ax[2].set_xlabel('(log-rate slope)$^2$'); ax[2].set_ylabel(r'$\ln\nu$')
    ax[2].set_title(r'curvature-free identity, slope $=-A/2b^2$', fontsize=8)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# --------------------------------------------------- out-of-sample delay slope
# Cell labels are duration/coherence ("0.3/c12" = 0.3 s duration, 12% coherence);
# the out-of-sample figures group by duration and spell coherence out instead, since
# that compact form reads poorly on an axis.
DUR_LAB = {0.3: '0.3 s', 2.0: '2.0 s'}
COH_LAB = {0.0: '0%', 0.12: '12%', 0.50: '50%'}
OOS_X = [0, 1, 2, 3.6, 4.6, 5.6]           # gap between the two duration clusters


def _oos_predictions(p1_behav, p2_behav, n=30000):
    """Simulate both behavior-only-calibrated models at every cell and return, per
    cell, both the full 20-bin delay-period rate profile and its log-rate slope --
    shared by fig_oos_slope and fig_oos_rate_profile so each model is only run once."""
    TB = np.arange(20) * 0.1 + 0.05; TB = TB - TB.mean()

    def run(p, model, seed0):
        profiles = np.zeros((6, 20)); slopes = np.zeros(6)
        for k, (c, D) in enumerate(CELLS):
            s = simulate(p, c, D, n, model, seed=seed0 + k)
            m = measures(s)
            pr = rate_profile(m, D + 1, D + 3, 20)
            profiles[k] = pr
            slopes[k] = poisson_slope(pr, TB)[1] if (pr > 0).sum() >= 4 else np.nan
        return profiles, slopes

    prof1, s1 = run(p1_behav, 1, 9000)
    prof2, s2 = run(p2_behav, 2, 9100)
    return dict(profile1=prof1, profile2=prof2, slope1=s1, slope2=s2)


def fig_oos_slope(oos, E, fname='plots/fig_oos_slope.pdf'):
    """Out-of-sample validation (Results, Simulation validation / Delay-locked): each
    model calibrated on behavior alone (results/calibq_m{1,2}.json), then simulated,
    without recalibration, at the delay-period partial-activation log-rate slope --
    a quantity neither calibration was allowed to see. Returns the numbers that fill
    the manuscript's SLOPE1/SLOPE2/ERR1/ERR2/ERRFAC placeholders."""
    s1, s2 = oos['slope1'], oos['slope2']
    obs_persub = E['slope']                              # (nsub, 6)
    obs_mean = np.nanmean(obs_persub, axis=0)
    obs_n = np.sum(np.isfinite(obs_persub), axis=0)
    obs_se = np.nanstd(obs_persub, axis=0, ddof=1) / np.sqrt(np.maximum(obs_n, 1))
    err1 = float(np.mean(np.abs(s1 - obs_mean) / obs_se))
    err2 = float(np.mean(np.abs(s2 - obs_mean) / obs_se))

    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    ax.errorbar(OOS_X, obs_mean, yerr=obs_se, fmt='o', color='k', ms=5, capsize=3,
                label='observed (±1 SE)', zorder=3)
    ax.plot(OOS_X, s1, 's--', color='#c1272d', ms=5, label='Model 1 (behavior-only fit)')
    ax.plot(OOS_X, s2, '^-', color='#3182bd', ms=5, label='Model 2 (behavior-only fit)')
    ax.axhline(0, color='0.7', lw=.7, zorder=0)
    ax.axvline(2.8, color='0.85', lw=.7, zorder=0)
    ax.set_xticks(OOS_X)
    ax.set_xticklabels([COH_LAB[c] for c, D in CELLS], fontsize=7)
    ax.set_xlabel('motion coherence, grouped by stimulus duration', fontsize=7.5)
    for xc, D in [(1, 0.3), (4.6, 2.0)]:
        ax.text(xc, ax.get_ylim()[1] * 0.97, f'duration = {DUR_LAB[D]}', ha='center',
                fontsize=7, color='0.35')
    ax.set_ylabel('delay-period log-rate slope (s$^{-1}$)')
    ax.set_title('Out-of-sample: neither calibration saw this quantity', fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1.0))
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)

    return dict(obs_mean=obs_mean.tolist(), obs_se=obs_se.tolist(),
                model1_slope=s1.tolist(), model2_slope=s2.tolist(),
                err1=err1, err2=err2, errfac=err1 / err2)


def fig_oos_rate_profile(oos, E, fname='plots/fig_oos_rate_profile.pdf'):
    """The delay-period partial-activation rate itself, bin by bin (not just its
    slope): observed vs. each behavior-only-calibrated model's out-of-sample
    prediction, one panel per cell, arranged duration (rows) x coherence (columns)."""
    tb = np.arange(20) * 0.1 + 1.05                      # time from stimulus offset (s)
    obs = E['dly'].mean(0)                                # (6, 20), across participants
    fig, ax = plt.subplots(2, 3, figsize=(8.4, 4.6), sharex=True, sharey=True)
    for k, (c, D) in enumerate(CELLS):
        row = 0 if D == 0.3 else 1
        col = [0.0, 0.12, 0.50].index(c)
        a = ax[row, col]
        a.plot(tb, obs[k], 'o-', color='k', ms=2.5, lw=.9, label='observed')
        a.plot(tb, oos['profile1'][k], 's--', color='#c1272d', ms=2.5, lw=.9,
               label='Model 1 (behavior-only)')
        a.plot(tb, oos['profile2'][k], '^-', color='#3182bd', ms=2.5, lw=.9,
               label='Model 2 (behavior-only)')
        a.set_title(f'duration {DUR_LAB[D]}, coherence {COH_LAB[c]}', fontsize=7.5)
    for a in ax[1]:
        a.set_xlabel('time from stimulus offset (s)', fontsize=7.5)
    for a in ax[:, 0]:
        a.set_ylabel('PMA rate per 0.1 s bin', fontsize=7.5)
    ax[0, 2].legend(fontsize=6.5, frameon=False, loc='upper left')
    fig.suptitle('Out-of-sample delay-period rate: neither calibration saw this curve',
                 fontsize=8.5, y=1.02)
    fig.savefig(fname, bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- Model 1
if __name__ == '__main__':
    P2 = json.load(open('results/fit_model2.json')); P1 = json.load(open('results/fit_model1.json'))
    E = build()
    fig_trajectories(P2, 2, fname='plots/fig_traj_m2.pdf')
    fig_trajectories(P1, 1, fname='plots/fig_traj_m1.pdf')
    fig_single(P2)
    fig_rates(P2, E)
    fig_shape(E)

    p1b = json.load(open('results/calibq_m1.json'))
    p2b = json.load(open('results/calibq_m2.json'))
    oos_raw = _oos_predictions(p1b, p2b)
    fig_oos_rate_profile(oos_raw, E)
    oos = fig_oos_slope(oos_raw, E)
    disc1 = json.load(open('results/ckptq_m1.json'))['obj']
    disc2 = json.load(open('results/ckptq_m2.json'))['obj']
    json.dump(dict(disc1=disc1, disc2=disc2, **oos),
              open('results/oos_slope_check.json', 'w'), indent=1)
    print('figures written')
    print(f"\nDISC1 (model 1 behavior-only G2) = {disc1:.2f}")
    print(f"DISC2 (model 2 behavior-only G2) = {disc2:.2f}")
    print(f"SLOPE1 (model 1 predicted slopes) = {oos['model1_slope']}")
    print(f"SLOPE2 (model 2 predicted slopes) = {oos['model2_slope']}")
    print(f"observed slope range = {min(oos['obs_mean']):.3f} to {max(oos['obs_mean']):.3f} s^-1")
    print(f"ERR1 = {oos['err1']:.2f}, ERR2 = {oos['err2']:.2f}, ERRFAC = {oos['errfac']:.2f}")
