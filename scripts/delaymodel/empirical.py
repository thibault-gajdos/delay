"""
empirical.py -- the trials of one participant, read from the OSF archive (86cbw, path in
DELAY_DATA): condition, accuracy, response time, the times of stimulus onset and offset,
of the response signal, of the response-related muscle activation (the last EMG onset in
the responding channel before the button press) and of the response, and the onsets of
the partial muscle activations (every other EMG onset in the responding channel before
the response-related activation).  Trials with no response, a response before the
response signal, or a response time below 0.2 s are excluded.
"""
import numpy as np, pandas as pd, os

ROOT = os.path.expanduser(os.environ.get('DELAY_DATA', '/home/claude/data'))
SEG = os.path.join(ROOT, 'data/emg_segmented')
BEH = os.path.join(ROOT, 'data/behavioral')
SF = 1024.0
TRIG_STIM = {21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36}
TRIG_RS = {121, 122, 123, 124, 125, 126, 131, 132, 133, 134, 135, 136}
COH = {0.0: 0, 0.12: 1, 0.5: 2}


def trials(sid):
    b = np.load(os.path.join(BEH, f's{sid}.npy'), allow_pickle=True)
    ev = pd.read_csv(os.path.join(SEG, f's{sid}_segmentedEvts.csv'))
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
        # response-related activation: last onset (350) in the responding channel before
        # the response
        pre = np.flatnonzero((code[:i_rp] == 350.0) & (chan[:i_rp] == rch))
        if pre.size == 0:
            continue
        i_burst = pre[-1]
        t = (smp - smp[i_on]) / SF                       # time from stimulus onset
        on_sel = np.flatnonzero((code == 350.0) & (chan == rch))
        pma = on_sel[on_sel < i_burst]                   # PMA onsets, responding channel
        out.append(dict(
            coh=COH[coh[i]], dur=0 if dur[i] == 18 else 1, acc=acc[i], rt=rt[i],
            t_off=t[i_off], t_rs=t[i_rs], t_resp=t[i_rp], t_burst=t[i_burst],
            pma_t=t[pma]))
    return out
