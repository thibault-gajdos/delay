"""fig_predictions_m2.py -- Model 2's predictions at the representative parameter vector.

The exploratory counterpart of the manuscript's Figure 1: instead of one mean trajectory
per admissible parameter set, ONE TRACE PER CONDITION at a single vector, so that the
condition ordering the model predicts can be read off directly.  Model 2 only: the two
architectures differ in beta alone, which contributes identically in all six cells, so
their predictions for behaviour, PMET and MET are the same.

Top row:    the gated motor-preparatory signal y(t), time-locked to the four trial events.
Bottom row: mean accuracy, response time, premotor and motor execution time by condition
            -- the simulated counterpart of the manuscript's Figure 3.

The chronometric panels are drawn on a common 60 ms span (40 ms for MET), the scale of the
observed duration effect, so that an effect the model produces at a tenth of the observed
size reads as small rather than as a magnified wiggle.

VECTOR selects the parameter vector:
  'central'   results/representative.json -- the admissible draw closest to the
              coordinate-wise median of the sampled parameters, with beta at the Model 2
              admissible median.  Central by construction, and not chosen to resemble the
              data: at this vector the behavioural effects are in the right direction but
              about a tenth of the observed size.
  'rounded'   ROUNDED (below) -- round numbers, chosen to sit near the data.  Selected by
              scoring every admissible draw against the six observed condition means for
              accuracy, RT, PMET, MET and delay-period rate, then rounding the best one to
              one or two significant figures and re-checking that the rounded vector is
              still admissible.  An illustration, not a fit: nothing in the paper's
              argument uses it.
  an integer  that draw index in the sweep dump, with its own beta.

Writes plots/fig_predictions_m2_<VECTOR>.{pdf,png}.
Usage: python3 fig_predictions_m2.py [NTRIALS] [VECTOR]
"""
import numpy as np, json, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from delaymodel.sim import simulate, measures, DT
from delaymodel.targets import CELLS

HERE = os.path.dirname(os.path.abspath(__file__))
NTR = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
VECTOR = sys.argv[2] if len(sys.argv) > 2 else 'central'
MODEL = 2

# colour: family by duration (0.3 s cyan, 2 s red), darkness by coherence
COL = {(0.0, 0.3): '#8fdfea', (0.12, 0.3): '#22a7bd', (0.5, 0.3): '#0d5f6e',
       (0.0, 2.0): '#f7b0b0', (0.12, 2.0): '#e24a4a', (0.5, 2.0): '#8f1717'}
DCOL = {0.3: '#22a7bd', 2.0: '#e24a4a'}
CC = [0.0, 0.12, 0.5]

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'lines.linewidth': 1.3})

ROUNDED = dict(sigx=0.1, Lmotion=0.2, LRS=0.2, delay=3.0, Lmotor=0.093,
               g=0.16, theta=0.06, kappa=12.0, SDe=0.05, k=0.9, lam=40.0,
               muRS=0.8, r=0.06, sigq=0.03)
ROUNDED_BETA = 0.04

if VECTOR == 'rounded':
    P = dict(ROUNDED, alpha=ROUNDED_BETA / 3.0)
    tag = 'rounded vector near the data'
elif VECTOR == 'central':
    rep = json.load(open(os.path.join(HERE, 'results/representative.json')))
    P = dict(rep['params'], alpha=rep['alpha'])
    tag = f"representative draw {rep['draw']}, beta at the Model 2 admissible median"
else:
    d = int(VECTOR)
    rec = next(json.loads(l) for l in open(os.path.join(HERE, 'results/sweep_screened.jsonl'))
               if json.loads(l)['draw'] == d)
    P = dict(rec['params'], alpha=rec['alpha'])
    tag = f"draw {d}, its own beta"
print(tag)

TG_ON = np.linspace(0.0, 2.0, 400)        # from stimulus onset
TG_OFF = np.linspace(0.0, 3.0, 400)       # from stimulus offset
TG_RS = np.linspace(0.0, 0.30, 120)       # from the response signal
L_RESP = int(0.20 / DT)                   # 200 ms before the response

traces, behav = {}, {}
for k, (c, D) in enumerate(CELLS):
    s_ = simulate(P, c, D, NTR, MODEL, seed=4400 + 17 * k)
    m = measures(s_)
    v = np.flatnonzero(m['valid'])
    Y = s_['Z'] + P['g']              # the signal itself, before the gate is subtracted
    ri = s_['resp_i']
    nT = Y.shape[1]
    idx = np.arange(nT)[None, :]      # censor each trial at its own response, then average
    mu = np.nanmean(np.where(idx <= ri[v][:, None], Y[v], np.nan), 0)
    t = np.arange(nT) * DT
    R = np.full((v.size, L_RESP), np.nan)
    for a_, j in enumerate(v):
        if ri[j] >= L_RESP:
            R[a_] = Y[j, ri[j] - L_RESP:ri[j]]
    traces[(c, D)] = dict(
        on=np.interp(TG_ON, t, mu, left=np.nan, right=np.nan),
        off=np.interp(TG_OFF, t - D, mu, left=np.nan, right=np.nan),
        rs=np.interp(TG_RS, t - (D + P['delay']), mu, left=np.nan, right=np.nan),
        resp=np.nanmean(R, 0))
    behav[(c, D)] = {q: float(np.nanmean(m[q][v])) for q in ('acc', 'rt', 'pmet', 'met')}
    del s_, Y, R
    b = behav[(c, D)]
    print(f"c={c:>5} D={D}  n_valid={v.size:5d}  acc={b['acc']:.3f} "
          f"rt={1000*b['rt']:6.1f} pmet={1000*b['pmet']:6.1f} met={1000*b['met']:5.1f} ms",
          flush=True)

# ----------------------------------------------------------------------- figure
fig = plt.figure(figsize=(8.6, 5.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.55, wspace=0.42)
gsA = gs[0, :].subgridspec(1, 4, width_ratios=[2.0, 3.0, 0.62, 0.62], wspace=0.16)
axA = [fig.add_subplot(gsA[0, i]) for i in range(4)]
axB = [fig.add_subplot(gs[1, i]) for i in range(4)]

grids = [TG_ON, TG_OFF, TG_RS, np.arange(-L_RESP, 0) * DT]
for (c, D) in CELLS:
    for a, g, key in zip(axA, grids, ['on', 'off', 'rs', 'resp']):
        a.plot(g, traces[(c, D)][key], color=COL[(c, D)])

# the delay windows and the two response windows get their own vertical scales, as in
# the empirical figure: the signal crosses the gate only after the response signal, and a
# single scale would flatten the delay ramp into the axis
for a in axA[:2]:
    a.set_ylim(0, P['g'] * 1.18)
    a.axhline(P['g'], color='0.35', ls=':', lw=.8)
for a in axA[2:]:
    a.set_ylim(0, (P['g'] + P['r']) * 1.05)
    a.axhline(P['g'], color='0.35', ls=':', lw=.8)
    a.axhline(P['g'] + P['r'], color='0.55', ls=':', lw=.8)
axA[0].axvline(0.3, color=DCOL[0.3], ls='--', lw=.8)
axA[0].axvline(2.0, color=DCOL[2.0], ls='--', lw=.8)
axA[0].set_ylabel('predicted mean motor\npreparation $y(t)$  (a.u.)')
axA[0].text(0.06, P['g'] * 1.02, 'gate $g$', fontsize=6.5, color='0.35')
axA[2].text(0.01, (P['g'] + P['r']) * 0.93, 'criterion $g+r$', fontsize=6.5, color='0.55')
for a, lab in zip(axA, ['from stimulus\nonset (s)', 'from stimulus\noffset (s)',
                        'from RS (s)', 'from resp. (s)']):
    a.set_xlabel(lab)
axA[2].set_xticks([0, 0.3]); axA[3].set_xticks([-0.2, 0])
axA[1].set_yticklabels([])
axA[2].set_yticklabels([])
axA[3].yaxis.tick_right()
axA[0].set_title('A', loc='left', fontweight='bold', fontsize=11, pad=6)

handles = [plt.Line2D([], [], color=COL[(c, D)],
                      label=f'{c:.0%}, {"0.3" if D == 0.3 else "2"} s')
           for D in (0.3, 2.0) for c in CC]
axA[1].legend(handles=handles, fontsize=6, frameon=False, ncol=3, loc='lower right',
              bbox_to_anchor=(1.0, -0.02), handlelength=1.3, columnspacing=1.0,
              labelspacing=0.2)

PANELS = [('acc', 'predicted mean accuracy', None),
          ('rt', 'predicted mean RT (s)', 0.060),
          ('pmet', 'predicted mean PMET (s)', 0.060),
          ('met', 'predicted mean MET (s)', 0.040)]
for a, (key, lab, span) in zip(axB, PANELS):
    vals = []
    for D in (0.3, 2.0):
        y = [behav[(c, D)][key] for c in CC]
        vals += y
        a.plot([0, 1, 2], y, '-o', color=DCOL[D], ms=4,
               label=('0.3 s' if D == 0.3 else '2 s'))
    a.set_xticks([0, 1, 2]); a.set_xticklabels(['0', '12', '50'])
    a.set_xlim(-0.25, 2.25)
    a.set_xlabel('coherence (%)'); a.set_ylabel(lab)
    if span:
        mid = 0.5 * (min(vals) + max(vals))
        a.set_ylim(mid - span / 2, mid + span / 2)
axB[0].set_ylim(0.4, 1.02)
axB[0].set_title('B', loc='left', fontweight='bold', fontsize=11, pad=6)
axB[0].legend(fontsize=6.5, frameon=False, loc='lower right', title='duration',
              title_fontsize=6.5)

fig.text(0.5, 0.005, f'Model 2, {tag}; {NTR} trials per condition', ha='center',
         fontsize=6.5, color='0.45')

os.makedirs(os.path.join(HERE, 'plots'), exist_ok=True)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, f'plots/fig_predictions_m2_{VECTOR}.{ext}'),
                bbox_inches='tight', **({'dpi': 200} if ext == 'png' else {}))
print(f'wrote plots/fig_predictions_m2_{VECTOR}.pdf and .png')
