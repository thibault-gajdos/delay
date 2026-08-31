"""figures_fit.py -- fit-quality figures for both fits and both models.

Quantile-probability plots in the Ratcliff idiom (Ratcliff & Tuerlinckx 2002; the format
used by Dendauw et al. 2024): for each condition the .1/.3/.5/.7/.9 quantiles of the
correct responses are plotted at x = p(correct) and those of the errors at x = 1 -
p(correct), so accuracy and the whole shape of both distributions appear together.
Quantiles of the same level are joined across conditions.

  fig_qp_rt.pdf        fit A: response time, observed vs both fitted models
  fig_qp_fitB.pdf      fit B: PMET and MET, observed vs both fitted models
  fig_pma_delay.pdf    fit B: delay-window activation onset shape and count distribution
  fig_pma_windows.pdf  fit B: P(>=1 activation) in each analysis window

Fit A's errors are pooled across participants and fit B's are pooled across participants
AND cells (see build_targets_b.py), so fit B's error quantiles are a single distribution
and are drawn once, not per condition.
"""
import numpy as np, json, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.behav_sim import simulate_behav, behav_rt_split
import calib_fitB as B

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8.5, 'legend.fontsize': 7,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
C_OBS, C_M1, C_M2 = 'k', '#c1272d', '#3182bd'
Q = [10, 30, 50, 70, 90]
TA = json.load(open('results/targets_q.json'))
TB = json.load(open('results/targets_b.json'))
PA = {t: json.load(open(f'results/calibq_{t}.json')) for t in ('m1', 'm2')}
PB = {t: json.load(open(f'results/fitB_{t}.json')) for t in ('m1', 'm2')}
LABS = [c['label'] for c in TA['cells']]


def fitA_pred(p, model, n=20000):
    """predicted p(correct), correct/error RT quantiles and premature rate, per cell"""
    out = []
    for k, cell in enumerate(TA['cells']):
        o = simulate_behav(p, cell['coh'], cell['dur'], n, model=model, seed=2000 + k)
        rc, ri, nprem = behav_rt_split(o)
        nv = max(rc.size + ri.size, 1)
        out.append(dict(p_cor=rc.size / nv,
                        q_cor=np.percentile(rc, Q) if rc.size > 20 else np.full(5, np.nan),
                        q_inc=np.percentile(ri, Q) if ri.size > 20 else np.full(5, np.nan),
                        p_prem=nprem / max(nv + nprem, 1)))
    return out


def fitB_pred(p, model, n=8000):
    q = dict(p); q.pop('alpha', None) if model == 1 else None
    out, pm_i, mt_i = [], [], []
    for k, cell in enumerate(TB['cells']):
        st = B.sim_stats(q, cell['coh'], cell['dur'], n, model, seed=5000 + 100 * k)
        nv = max(st['n_cor'] + st['n_inc'], 1)
        out.append(dict(p_cor=st['n_cor'] / nv,
                        q_pmet=np.percentile(st['pmet_c'], Q) if st['pmet_c'].size > 20
                               else np.full(5, np.nan),
                        q_met=np.percentile(st['met_c'], Q) if st['met_c'].size > 20
                              else np.full(5, np.nan),
                        count=st['count_p'], shape=st['shape_p'], p_any=st['p_any']))
        pm_i.append(st['pmet_i']); mt_i.append(st['met_i'])
    pooled = dict(pmet=np.percentile(np.concatenate(pm_i), Q),
                  met=np.percentile(np.concatenate(mt_i), Q))
    return out, pooled


def _srt(x, qs, i):
    """the i-th quantile across conditions, ordered by response proportion"""
    pair = sorted((xx, q[i]) for xx, q in zip(x, qs) if np.isfinite(q[i]))
    return [p[0] for p in pair], [p[1] for p in pair]


def qp_panel(ax, x_obs, q_obs_cor, q_obs_inc, preds, ylab, title, xlim=(-0.02, 1.02)):
    """one quantile-probability panel.  preds: list of (label, colour, marker, x, qcor, qinc)"""
    for i in range(5):
        ax.plot(*_srt(x_obs, q_obs_cor, i), 'o-', color=C_OBS, ms=4, lw=1.0,
                zorder=3, label='observed' if i == 0 else None)
        if q_obs_inc is not None:
            ax.plot(*_srt([1 - x for x in x_obs], q_obs_inc, i), 'o-', color=C_OBS,
                    ms=4, lw=1.0, mfc='w', zorder=3)
    for lab, col, mk, xs, qc, qi in preds:
        for i in range(5):
            ax.plot(*_srt(xs, qc, i), mk + '--', color=col, ms=3.5, lw=.8, alpha=.8,
                    label=lab if i == 0 else None)
            if qi is not None:
                ax.plot(*_srt([1 - x for x in xs], qi, i), mk + '--', color=col,
                        ms=3.5, lw=.8, alpha=.8, mfc='w')
    ax.set_xlabel('response proportion'); ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=8.5)
    ax.set_xlim(*xlim)


def _srt(x, qs, i):
    """the i-th quantile across conditions, ordered by response proportion"""
    pair = sorted((xx, q[i]) for xx, q in zip(x, qs) if np.isfinite(q[i]))
    return [p[0] for p in pair], [p[1] for p in pair]


def qp_grid(fname, title, obs_x, obs_qc, obs_qi, pred, ylab, pooled=None):
    """Quantile-probability, one ROW per model and one COLUMN per stimulus duration.

    Splitting this way is what makes it readable: the full six-condition version puts
    3 series x 5 quantiles x 2 response types = 30 lines in one frame.  Here each panel
    carries one model against the data over three coherences, and correct and error
    responses fall on opposite halves of the axis.

    obs_x/obs_qc/obs_qi are per-cell (6).  pred is {model_tag: (x, qcor, qinc)}.
    `pooled`, if given, is (obs5, {tag: pred5}) for errors pooled over all conditions.
    """
    ncol = 2 + (1 if pooled else 0)
    fig, axs = plt.subplots(2, ncol, figsize=(4.1 * ncol, 6.4),
                            gridspec_kw=dict(width_ratios=[1, 1] + ([0.42] if pooled else [])))
    for r, (tag, col, mk, name) in enumerate([('m1', C_M1, 's', 'Model 1'),
                                              ('m2', C_M2, '^', 'Model 2')]):
        px, pqc, pqi = pred[tag]
        for cix, (D, dlab) in enumerate([(0.3, 'stimulus 0.3 s'), (2.0, 'stimulus 2.0 s')]):
            ax = axs[r, cix]
            sel = [k for k, c in enumerate(TA['cells']) if c['dur'] == D]
            if obs_qi is not None:
                ax.axvspan(-0.02, 0.5, color='0.94', zorder=0)
                ax.axvline(0.5, color='0.75', lw=.7, zorder=1)
            ox = [obs_x[k] for k in sel]
            for i in range(5):
                ax.plot(*_srt(ox, [obs_qc[k] for k in sel], i), 'o-', color=C_OBS,
                        ms=4.5, lw=1.1, zorder=4, label='observed' if i == 0 else None)
                ax.plot(*_srt(px[sel[0]:sel[-1] + 1] if False else [px[k] for k in sel],
                              [pqc[k] for k in sel], i), mk + '--', color=col, ms=4,
                        lw=.9, zorder=3, label=name if i == 0 else None)
                if obs_qi is not None:
                    ax.plot(*_srt([1 - obs_x[k] for k in sel], [obs_qi[k] for k in sel], i),
                            'o-', color=C_OBS, ms=4.5, lw=1.1, mfc='w', zorder=4)
                    ax.plot(*_srt([1 - px[k] for k in sel], [pqi[k] for k in sel], i),
                            mk + '--', color=col, ms=4, lw=.9, mfc='w', zorder=3)
            # fit B has no per-condition error quantiles (they are pooled over all
            # conditions and drawn in their own panel), so the error half would be blank
            ax.set_xlim(-0.02 if obs_qi is not None else 0.45, 1.02)
            if r == 0:
                ax.set_title(dlab, fontsize=9)
            if r == 1:
                ax.set_xlabel('response proportion')
            if cix == 0:
                ax.set_ylabel(f'{name}\n{ylab}')
                if r == 0:
                    ax.legend(frameon=False, loc='upper center', ncol=2)
            if obs_qi is not None:
                ax.text(0.25, 0.965, 'errors', transform=ax.transAxes, ha='center',
                        fontsize=7, color='0.45')
                ax.text(0.75, 0.965, 'correct', transform=ax.transAxes, ha='center',
                        fontsize=7, color='0.45')
            else:
                ax.text(0.5, 0.965, 'correct responses', transform=ax.transAxes,
                        ha='center', fontsize=7, color='0.45')
        if pooled:
            ax = axs[r, 2]
            o5, pr = pooled
            ax.plot([0]*5, o5, 'o', color=C_OBS, ms=5, label='observed')
            ax.plot([1]*5, pr[tag], mk, color=col, ms=5, label=name)
            for i in range(5):
                ax.plot([0, 1], [o5[i], pr[tag][i]], '-', color='0.8', lw=.7, zorder=0)
            ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1])
            ax.set_xticklabels(['obs', name.split()[-1]], fontsize=7.5)
            ax.set_title('errors,\npooled over\nall conditions', fontsize=7.5)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f'plots/{fname}', bbox_inches='tight'); plt.close(fig)


# ---------------------------------------------------------------- 1. RT, fit A
obsA = TA['cells']
p1, p2 = fitA_pred(PA['m1'], 1), fitA_pred(PA['m2'], 2)
qp_grid('fig_qp_rt.pdf', 'Fit A (behaviour only): response-time quantile-probability',
        [c['p_correct'] for c in obsA],
        [c['correct_quantiles'] for c in obsA],
        [c['incorrect_quantiles'] for c in obsA],
        {'m1': ([d['p_cor'] for d in p1], [d['q_cor'] for d in p1], [d['q_inc'] for d in p1]),
         'm2': ([d['p_cor'] for d in p2], [d['q_cor'] for d in p2], [d['q_inc'] for d in p2])},
        'response time (s)')

# ---------------------------------------------------------------- 2. PMET / MET, fit B
b1, pool1 = fitB_pred(PB['m1'], 1)
b2, pool2 = fitB_pred(PB['m2'], 2)
obsB = TB['cells']
for key, okey, pkey, fname, lab in [
        ('q_pmet', 'pmet_quantiles_correct', 'pmet', 'fig_qp_pmet.pdf', 'premotor time PMET (s)'),
        ('q_met', 'met_quantiles_correct', 'met', 'fig_qp_met.pdf', 'motor execution time MET (s)')]:
    qp_grid(fname, f'Fit B (all data): {lab}',
            [c['p_correct'] for c in obsB], [c[okey] for c in obsB], None,
            {'m1': ([d['p_cor'] for d in b1], [d[key] for d in b1], None),
             'm2': ([d['p_cor'] for d in b2], [d[key] for d in b2], None)},
            lab,
            pooled=(TB[f'{pkey}_quantiles_incorrect_pooled'],
                    {'m1': pool1[pkey], 'm2': pool2[pkey]}))

# ---------------------------------------------------------------- 3. delay PMA
tb = np.arange(TB['n_shape_bins']) * 0.4 + 1.2      # bin centres, s from stimulus offset
fig, axs = plt.subplots(2, 6, figsize=(12.4, 4.4), sharex='row', sharey='row')
for k, cell in enumerate(obsB):
    a = axs[0, k]
    a.plot(tb, cell['delay_shape_p'], 'o-', color=C_OBS, ms=3.5, lw=1, label='observed')
    a.plot(tb, b1[k]['shape'], 's--', color=C_M1, ms=3.5, lw=.9, label='Model 1')
    a.plot(tb, b2[k]['shape'], '^--', color=C_M2, ms=3.5, lw=.9, label='Model 2')
    a.axhline(1.0 / TB['n_shape_bins'], color='0.8', lw=.7, zorder=0)
    a.set_title(cell['label'], fontsize=7.5)
    if k == 0:
        a.set_ylabel('share of activations')
    b = axs[1, k]
    xc = np.arange(4)
    b.plot(xc, cell['delay_count_p'], 'o-', color=C_OBS, ms=3.5, lw=1)
    b.plot(xc, b1[k]['count'], 's--', color=C_M1, ms=3.5, lw=.9)
    b.plot(xc, b2[k]['count'], '^--', color=C_M2, ms=3.5, lw=.9)
    b.set_xticks(xc); b.set_xticklabels(['0', '1', '2', '3+'])
    b.set_xlabel('activations / trial')
    if k == 0:
        b.set_ylabel('proportion of trials')
axs[0, 0].legend(frameon=False, fontsize=6.5)
for a in axs[0, :]:
    a.set_xlabel('time from stimulus offset (s)', fontsize=6.8)
fig.subplots_adjust(hspace=0.62)
fig.suptitle('Fit B: delay-window partial activations -- onset shape (top) and count per '
             'trial (bottom)', fontsize=9)
fig.savefig('plots/fig_pma_delay.pdf', bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- 4. P(>=1) per window
fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.2), sharex=True)
x = np.arange(6)
for j, (key, name) in enumerate([('p_any_stim', 'stimulus window'),
                                 ('p_any_rs', 'RS window'),
                                 ('p_any_resp', 'response window')]):
    a = axs[j]
    a.plot(x, [c[key] for c in obsB], 'o-', color=C_OBS, ms=4, lw=1, label='observed')
    a.plot(x, [d['p_any'][j] for d in b1], 's--', color=C_M1, ms=4, lw=.9, label='Model 1')
    a.plot(x, [d['p_any'][j] for d in b2], '^--', color=C_M2, ms=4, lw=.9, label='Model 2')
    a.set_xticks(x); a.set_xticklabels(LABS, rotation=45, ha='right', fontsize=6.5)
    a.set_title(name, fontsize=8.5)
    if j == 0:
        a.set_ylabel('P(at least one partial activation)'); a.legend(frameon=False)
fig.savefig('plots/fig_pma_windows.pdf', bbox_inches='tight'); plt.close(fig)

print("written: fig_qp_rt.pdf, fig_qp_pmet.pdf, fig_qp_met.pdf, fig_pma_delay.pdf, fig_pma_windows.pdf")
