"""
empirical.py -- every empirical quantity used in the note, extracted from the OSF
archive (86cbw).  Counting follows analysis_PMA/compute_PMA.py: a partial activation
is an EMG onset (code 350) in the channel of the responding effector, occurring
before the onset of the response burst.

Windows (manuscript p. 25, lines 516-521):
    stimulus  0 to 2 s after stimulus onset          20 bins
    delay     1 to 3 s after stimulus offset         20 bins
    RS        0 to 0.2 s after RS onset               2 bins
    response  -0.2 to 0 s before the button press     2 bins
"""
import numpy as np, pandas as pd, os

ROOT = os.path.expanduser(os.environ.get('DELAY_DATA', '/home/claude/data'))
SEG = os.path.join(ROOT, 'data/emg_segmented')
BEH = os.path.join(ROOT, 'data/behavioral')
SF = 1024.0
TRIG_STIM = {21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36}
TRIG_RS = {121, 122, 123, 124, 125, 126, 131, 132, 133, 134, 135, 136}
COH = {0.0: 0, 0.12: 1, 0.5: 2}


def load_subject(sid):
    b = np.load(os.path.join(BEH, f's{sid}.npy'), allow_pickle=True)
    ev = pd.read_csv(os.path.join(SEG, f's{sid}_segmentedEvts.csv'))
    return b, ev


def trials(sid):
    b, ev = load_subject(sid)
    coh = b[:, 5].astype(float); dur = b[:, 6].astype(float).astype(int)
    rmom = b[:, 10].astype(float).astype(int); acc = b[:, 11].astype(float).astype(int)
    rt = b[:, 12].astype(float); trg = b[:, 14].astype(float).astype(int)
    out = []
    for idx, gp in ev.groupby('trialIdx'):
        i = int(idx)
        if i >= len(rt):
            continue
        if rmom[i] != 1 or acc[i] == -1 or rt[i] < 0.2:
            continue
        code = gp['code'].to_numpy(); smp = gp['sample'].to_numpy(float)
        chan = gp['chan'].to_numpy()
        if 99.0 in code:
            continue
        try:
            i_on = np.flatnonzero(np.isin(code, list(TRIG_STIM)))[0]
            i_off = np.flatnonzero(code == 1.0)[0]
            i_rs = np.flatnonzero(np.isin(code, list(TRIG_RS)))[0]
            i_rp = np.flatnonzero(code == float(trg[i]))[0]
        except IndexError:
            continue
        rch = 0 if trg[i] in (100, 101) else 1
        # response burst: last onset (350) in the responding channel before the response
        pre = np.flatnonzero((code[:i_rp] == 350.0) & (chan[:i_rp] == rch))
        if pre.size == 0:
            continue
        i_burst = pre[-1]
        t = (smp - smp[i_on]) / SF                       # time from stimulus onset
        on_sel = np.flatnonzero((code == 350.0) & (chan == rch))
        off_sel = np.flatnonzero((code == 351.0) & (chan == rch))
        pma = on_sel[on_sel < i_burst]                   # PMA onsets, selected channel
        on_opp = np.flatnonzero((code == 350.0) & (chan == 1 - rch))
        # duration of each PMA: matching 351 in the same channel after the 350
        dur_pma = []
        for k in pma:
            nxt = off_sel[off_sel > k]
            if nxt.size:
                dur_pma.append(t[nxt[0]] - t[k])
        out.append(dict(
            coh=COH[coh[i]], dur=0 if dur[i] == 18 else 1, acc=acc[i], rt=rt[i],
            t_off=t[i_off], t_rs=t[i_rs], t_resp=t[i_rp], t_burst=t[i_burst],
            pma_t=t[pma], pma_d=np.array(dur_pma),
            opp_t=t[on_opp[on_opp < i_burst]],
            pma_off=np.array([t[off_sel[off_sel > k][0]] for k in pma
                              if off_sel[off_sel > k].size]),
            burst_off=(t[off_sel[off_sel > i_burst][0]]
                       if off_sel[off_sel > i_burst].size else np.nan),
        ))
    return out


def build(nsub=22):
    stim = np.zeros((nsub, 2, 3, 20)); delay = np.zeros((nsub, 2, 3, 20))
    rs = np.zeros((nsub, 2, 3, 2)); resp = np.zeros((nsub, 2, 3, 2))
    ntr = np.zeros((nsub, 2, 3)); opp = np.zeros((nsub, 2, 3, 20))
    beh = {k: np.full((nsub, 2, 3), np.nan) for k in ['acc', 'rt', 'pmet', 'met']}
    durs, gaps, tail = [], [], []
    for s in range(nsub):
        T = trials(s + 1)
        acc = {}; rt = {}; pmet = {}; met = {}
        for tr in T:
            d, c = tr['dur'], tr['coh']
            ntr[s, d, c] += 1
            # stimulus period, 0-2 s from onset, truncated at the response burst
            v = tr['pma_t']
            k = v[(v >= 0) & (v < 2.0) & (v < tr['t_burst'])]
            stim[s, d, c] += np.histogram(k, bins=20, range=(0, 2))[0]
            o = tr['opp_t']; o = o[(o >= 0) & (o < 2.0)]
            opp[s, d, c] += np.histogram(o, bins=20, range=(0, 2))[0]
            # delay, 1-3 s from offset
            k = v[(v >= tr['t_off'] + 1) & (v < tr['t_off'] + 3)] - tr['t_off'] - 1
            delay[s, d, c] += np.histogram(k, bins=20, range=(0, 2))[0]
            # RS, 0-0.2 s
            k = v[(v >= tr['t_rs']) & (v < tr['t_rs'] + 0.2)] - tr['t_rs']
            rs[s, d, c] += np.histogram(k, bins=2, range=(0, 0.2))[0]
            # response, -0.2 to 0 before the press
            k = v[(v >= tr['t_resp'] - 0.2) & (v < tr['t_resp'])] - tr['t_resp'] + 0.2
            resp[s, d, c] += np.histogram(k, bins=2, range=(0, 0.2))[0]
            acc.setdefault((d, c), []).append(tr['acc'])
            rt.setdefault((d, c), []).append(tr['rt'])
            pmet.setdefault((d, c), []).append(tr['t_burst'] - tr['t_rs'])
            met.setdefault((d, c), []).append(tr['t_resp'] - tr['t_burst'])
            # burst durations in the delay, and gaps between consecutive activations
            m = (tr['pma_t'] >= tr['t_off'] + 1) & (tr['pma_t'] < tr['t_off'] + 3)
            if tr['pma_d'].size == tr['pma_t'].size:
                durs.extend(tr['pma_d'][m])
            if tr['pma_off'].size == tr['pma_t'].size and tr['pma_t'].size > 1:
                gaps.extend(tr['pma_t'][1:] - tr['pma_off'][:-1])
            # gap between the last PMA and the response burst
            if tr['pma_off'].size and tr['pma_t'].size == tr['pma_off'].size:
                tail.append(tr['t_burst'] - tr['pma_off'][-1])
        for (d, c), v in acc.items():
            beh['acc'][s, d, c] = np.mean(v); beh['rt'][s, d, c] = np.mean(rt[(d, c)])
            beh['pmet'][s, d, c] = np.mean(pmet[(d, c)])
            beh['met'][s, d, c] = np.mean(met[(d, c)])
    n = ntr[..., None]
    np.savez('empirical.npz', stim=stim / n, delay=delay / n, rs=rs / n,
             resp=resp / n, opp=opp / n, ntr=ntr,
             durs=np.array(durs), gaps=np.array(gaps), tail=np.array(tail), **beh)
    return dict(stim=stim / n, delay=delay / n, rs=rs / n, resp=resp / n,
                ntr=ntr, durs=np.array(durs), gaps=np.array(gaps),
                tail=np.array(tail), **beh)


if __name__ == '__main__':
    E = build()
    L = ['0.3 c00', '0.3 c12', '0.3 c50', '2.0 c00', '2.0 c12', '2.0 c50']
    print(f"valid trials: {int(E['ntr'].sum())}")
    print("\ncell        acc     rt    pmet    met   |  stim   delay    rs   resp  (PMA/0.1s bin)")
    for k, (d, c) in enumerate([(a, b) for a in range(2) for b in range(3)]):
        print(f"{L[k]:8s} {np.nanmean(E['acc'][:,d,c]):6.3f} "
              f"{np.nanmean(E['rt'][:,d,c]):6.3f} {np.nanmean(E['pmet'][:,d,c]):6.3f} "
              f"{np.nanmean(E['met'][:,d,c]):6.3f}   | "
              f"{E['stim'][:,d,c].mean():6.4f} {E['delay'][:,d,c].mean():6.4f} "
              f"{E['rs'][:,d,c].mean():6.4f} {E['resp'][:,d,c].mean():6.4f}")
    print(f"\npooled delay rate  : {E['delay'].mean():.4f} per 0.1 s bin")
    print(f"pooled stimulus    : {E['stim'].mean():.4f}")
    print(f"delay burst duration : n={E['durs'].size}, median {1000*np.median(E['durs']):.0f} ms, "
          f"mean {1000*E['durs'].mean():.0f} ms")
    g = E['gaps']; g = g[np.isfinite(g)]
    print(f"inter-activation gaps: n={g.size}, "
          f"< 25 ms: {100*(g<0.025).mean():.2f} %, < 50 ms: {100*(g<0.05).mean():.2f} %, "
          f"median {1000*np.median(g):.0f} ms")
    t = E['tail']; t = t[np.isfinite(t)]
    print(f"gap last PMA -> response burst: median {1000*np.median(t):.0f} ms, "
          f"< 50 ms: {100*(t<0.05).mean():.1f} %, < 25 ms: {100*(t<0.025).mean():.1f} %")
