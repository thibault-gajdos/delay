import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
rng = np.random.default_rng(20260826)
D = '/mnt/user-data/uploads/delay/data/86cbw-osfstorage-archive/analysis_behavioral_emg/'
COH = [0,12,50]; DUR = [300,2000]
COL = {300:'#1FCBD1', 2000:'#DE5859'}
sns.set_theme(style='darkgrid')
plt.rcParams.update({'font.size':11,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11})

def load(name, pref):
    df = pd.read_csv(D+name)
    A = np.empty((len(df), 2, 3))
    for j,d in enumerate(DUR):
        for k,c in enumerate(COH):
            A[:,j,k] = df[f'{pref}_coh{c}_{d}'].values
    return A

def morey_norm(A):
    """Cousineau-Morey normalisation: remove each participant's overall mean, add grand mean."""
    flat = A.reshape(A.shape[0], -1)
    return (flat - flat.mean(axis=1, keepdims=True) + flat.mean()).reshape(A.shape)

def morey_sem(A):
    n, M = A.shape[0], A.shape[1]*A.shape[2]
    N = morey_norm(A)
    return N.std(axis=0, ddof=1)/np.sqrt(n) * np.sqrt(M/(M-1))

def morey_boot(A, B=10000):
    """95% percentile bootstrap of the Cousineau-Morey normalised values, Morey-scaled about the mean."""
    n, M = A.shape[0], A.shape[1]*A.shape[2]
    N = morey_norm(A); mean = A.mean(axis=0)
    idx = rng.integers(0, n, size=(B, n))
    bs = N[idx].mean(axis=1)                       # B x 2 x 3
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    f = np.sqrt(M/(M-1)); nm = N.mean(axis=0)
    return mean - (nm-lo)*f, mean + (hi-nm)*f

panels = [('accuracy','accuracy','p(correct)','A'),
          ('mean_RT','mean_RT','Mean RT (s)','B'),
          ('mean_PMET','mean_PMT','Mean PMET (s)','C'),
          ('mean_MET','mean_MT','Mean MET (s)','D')]

fig, axes = plt.subplots(1, 4, figsize=(16.6, 4.138))
x = np.arange(3)
for ax,(f,pref,ylab,lab) in zip(axes, panels):
    A = load(f+'.csv', pref)
    for j,d in enumerate(DUR):
        for s in range(A.shape[0]):
            ax.plot(x, A[s,j,:], color=COL[d], alpha=0.22, lw=0.7, marker='.', ms=2.5, zorder=1)
    if f == 'accuracy':
        lo, hi = morey_boot(A)
    else:
        se = morey_sem(A); lo, hi = A.mean(axis=0)-se, A.mean(axis=0)+se
    for j,d in enumerate(DUR):
        m = A[:,j,:].mean(axis=0)
        ax.errorbar(x, m, yerr=[m-lo[j], hi[j]-m], color=COL[d], lw=2.6,
                    marker='o', ms=7, capsize=0, elinewidth=2.0, zorder=3,
                    label=f'{d} ms' if lab=='A' else None)
    ax.set_xticks(x); ax.set_xticklabels(['0','12','50'])
    ax.set_xlim(-0.25, 2.25); ax.set_xlabel('Coherence (%)'); ax.set_ylabel(ylab)
    ax.set_title(lab, loc='left', fontweight='bold', fontsize=14)
h,l = axes[0].get_legend_handles_labels()
fig.legend(h, l, title='Stimulus Duration', loc='center left', bbox_to_anchor=(0.895,0.5), frameon=False)
fig.tight_layout(rect=[0,0,0.885,1])
fig.savefig('fig3_new.png', dpi=200)
fig.savefig('fig3_new.tiff', dpi=370, pil_kwargs={'compression':'tiff_lzw'})
from PIL import Image
print('saved', Image.open('fig3_new.tiff').size)
