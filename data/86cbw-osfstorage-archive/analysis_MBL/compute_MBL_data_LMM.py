"""
MBL (8-30Hz) data 
"""

import numpy as np
import os
import re
import copy
import sys
import matplotlib.pyplot as plt
import pandas as pd
import arviz as az
import scipy.stats as st
from tqdm import tqdm          
np.set_printoptions(5, suppress=True)
az.style.use('arviz-whitegrid')
from matplotlib.ticker import FuncFormatter
current = os.path.join(os.getcwd()) 
toolbox_EMG_folder = r"C:\Users\erc\iCloudDrive\ACHARD_these\EMG_method\toolbox_EMG"
sys.path.extend([toolbox_EMG_folder]) # to add toolbox path to environment paths
import emgTools as emg
import evtTools as evt
from Events import Events
import mne



"""
Data behavior

|  | Columns Name           | Signification      | Var   | Details                                                                                                        |
|  |                        |                    |       |                                                                                                                |
|0 | subject_id             | ID Subject         | Int   | Subject number (22 valids)                                                                                     |
|1 | trial_tot              | Trial number tot   | Int   | Number of trial in total for all block (1 -> 960 trials) 150 trials/conditions )                               |
|2 | block                  | Block number       | Int   | Number of trial in total for all block.                                                                        |
|3 | block_type             | Type of block      | Str   | There is 300 ms or 2s.                                                                                         |
|4 | trial_block            | Trial number block | Int   | Number of trial for block (1 -> 96 trials)                                                                     |
|5 | condition_coherence    | Motion coherence   | Int   | Motion coherence of trial. Condition of coherence possible: 0, .12 or .50                                      |
|6 | stim_pres              | Time of pres. stim | Float | Time of stimulus duration during stimulus period: 300 ms or 2000ms                                             |
|7 | direction_motion       | Dots direction     | Int   | Direction of dots in RDM for the trial. Left is 180, Right is 0.                                               |
|8 | key_required           | Expected key       | Int   | The key response expected for a correct response (1-left or 2-right).                                          |
|9 | key_pressed            | Pressedey          | Int   | The key response pressed (1-left or 2-right).                                                                  |
|10| response_moment        | Time of response   | Int   | The moment when the participant pressed a key can be too early (0), at the right time (1), or to late (-1).    |
|11| accuracy               | Accuracy           | Int   | The accuracy of each subject for the task 1 if correct, 0 if incorrect, else -1.                               |
|12| rt                     | RT                 | Float | Response Time of participant (in ms) from RS                                                                   |
|13| trigg_stim             | Stimulus trigger   | Int   | The stimulus trigger indicate conditions select of trial (delay, coherence, direction of dots). See * below.   |
|14| trigg_resp             | Response trigger   | Int   | The response trigger indicate response subject (correct or not, left or right). See ** below.                  |



* Stimulus Onset Trigger for 3s delay
Trigger Stim code:

|                 | Stim Pres |  coherence  0.  |  coherence .12  |  coherence .50  |  
|                 |           |                 |                 |                 |  
| RS 300          |  300 ms   |   L 21 / R 24   |   L 22 / R 25   |   L 23 / R 26   |  
| RS 2000         |  2000 ms  |   L 31 / R 34   |   L 32 / R 35   |   L 33 / R 36   |  

Note. Response Signal trigger (after disapear of cross = trigg_stim + 100)



**** Trigger Response ****

| Response Moment | Trigger response | Details                      |
|                 |                  |                              |
| Too late  -1    | 250              | Too late, no key response    | 
|                 |                  |                              |
| Too early  0    | 200              | Incorrect left response      |
| Too early  0    | 202              | Incorrect right response     |
| Too early  0    | 201              | Correct left response        |
| Too early  0    | 203              | Correct right response       |
|                 |                  |                              |
| At time    1    | 100              | Incorrect left response      |
| At time    1    | 102              | Incorrect right response     |
| At time    1    | 101              | Correct left response        |
| At time    1    | 103              | Correct right response       |

"""


def index_nearest(array,value):#returns the index of the nearest value
    idx = (np.abs(array-value)).argmin()
    return idx


# ------------------------ Path ------------------------
path_behav = os.path.abspath(os.path.join(current, "..", "..", "data", "behavioral"))
path_tt = os.path.abspath(os.path.join(current, "..", "..", "data", "vector_trial_type"))
path_eeg = os.path.abspath(os.path.join(current, "..", "..", "processing_EEG", "data_saved"))
path_data_export = os.path.abspath(os.path.join(current, "..", "..", "data", "data_export"))

#  Load time_step ---
time_steps_onset  = np.load(os.path.join(path_eeg,"time_TF_onset_300ms.npy")) # remove padding, begin at zero
time_steps_offset = np.load(os.path.join(path_eeg,"time_TF_offset_300ms.npy"))
time_steps_rs     = np.load(os.path.join(path_eeg,"time_TF_rs_300ms.npy"))
time_steps_resp   = np.load(os.path.join(path_eeg,"time_TF_resp_300ms.npy"))

freq_stft = np.load(os.path.join(path_eeg,"freq_stft_300ms.npy"))
df_export = np.load(path_data_export + "\data_export.npy")


#  Create data containers for MBL voltage for long format ----
nb_subjects = 22
stim_duration = 2
coh = 3
eeg = 64
freqs = 6 # win 300ms
press = 3
    # 0 left press only (correct or incorrect)
    # 1 right left press only (correct or incorrect)
    # 2 all press take in account
    
possible_stim = [0.3, 2]
possible_coh = [0, 0.12, 0.50]

# Data for average by condition -----------------------------
# Electrods mapping was reorganised (all press as C3 electrod in controlateral response - like a right press which C3 is controlateral to response)
MB_onset_conditions = np.zeros((nb_subjects, stim_duration, coh, eeg,  freqs, len(time_steps_onset), press)) # (22, 2, 3, 64, 9831, 3)
MB_offset_conditions = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_offset), press))
MB_rs_conditions = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_rs),  press))
MB_resp_conditions = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_resp), press))

# No electrods mapping was reorganised (C3 is left H and C4 is right H)
MB_onset_conditions_no_reorg = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_onset),  press)) # (22, 2, 3, 64, 9831, 3 press, 13 freqs)
MB_offset_conditions_no_reorg = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_offset),  press))
MB_rs_conditions_no_reorg = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_rs), press))
MB_resp_conditions_no_reorg = np.zeros((nb_subjects, stim_duration, coh, eeg, freqs, len(time_steps_resp), press))




# Data for long format by condition -----------------------------
# Create Bin length (in sample) --- 
# e.g., 50ms is 0.05s * 1024Hz = 51.2 samples
# e.g., 100ms is 0.1s * 1024Hz = 102.4 samples

convert = 1 # to convert s in sample (sampling rate), e.g., 1024 Hz, else 1 if in s
end = 0.01 # to take account the end of the bin (+1 if sample, else 0.001 for different of 1ms)
len_bin = 0.1 * convert

# Create bin (in s to correspond epochTime array)---
max_after_RS = 0.5
bin_all = np.arange(-0.1, (2+3+ max_after_RS) * convert , len_bin)  # from stimulus onset 0 to +0.2s after RS onset
bin_stim = np.arange(-0.1, (2) * convert + end, len_bin)            # max 2s stimulus duration
bin_delay = np.arange(-0.1, (3) * convert + end , len_bin)          # 3s delay
bin_rs = np.arange(-0.1, (max_after_RS) * convert + end, len_bin)   # 400ms after RS
bin_resp = np.arange(-0.3, 0.3 * convert + end, len_bin)            # -0.3s before RS and +0.3s after RS


# Create Stimulus period data container
nb_columns = 9
n_total = 22 * 960 * len(bin_stim) 
arr_stim = np.zeros((n_total, nb_columns))
row_idx_stim = 0


# Create Stimulus period data container
n_total = 22 * 960 * len(bin_delay) 
arr_delay= np.zeros((n_total, nb_columns))
row_idx_delay = 0

# Create Stimulus period data container
n_total = 22 * 960 * len(bin_rs) 
arr_rs = np.zeros((n_total, nb_columns))
row_idx_rs = 0

# Create Stimulus period data container
n_total = 22 * 960 * len(bin_resp) 
arr_resp = np.zeros((n_total, nb_columns))
row_idx_resp = 0







# ============================================ Flip channel egg (all response in same hemisphere) ============================================
# --- For correct + incorrect ---
# if R resp = ipsi
# if L resp = ipsi


# ---- ch_names without EOG ----
sf = 1024
ch_names = np.array(['Fp1', 'AF7', 'AF3', 'F1', 'F3', 'F5', 'F7', 'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7', 'TP7',
            'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', 'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz',
            'Fpz', 'Fp2', 'AF8', 'AF4', 'AFz', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT8', 'FC6', 'FC4', 'FC2', 'FCz', 'Cz',
            'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO8', 'PO4', 'O2'])

ch_type = ['eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
           'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
           'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
           'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg', 'eeg',
           'eeg', 'eeg', 'eeg', 'eeg'] 

# ---- DataFrame with flipped colums and original eeg channel
# Dictionnaire de correspondance gauche-droite (flipping)
flipping_dict = {
    'Fp1': 'Fp2', 'Fp2': 'Fp1', 'AF7': 'AF8', 'AF8': 'AF7', 'AF3': 'AF4', 'AF4': 'AF3', 'F1': 'F2', 'F2': 'F1',
    'F3': 'F4', 'F4': 'F3', 'F5': 'F6', 'F6': 'F5', 'F7': 'F8', 'F8': 'F7', 'FT7': 'FT8', 'FT8': 'FT7',
    'FC5': 'FC6', 'FC6': 'FC5', 'FC3': 'FC4', 'FC4': 'FC3', 'FC1': 'FC2', 'FC2': 'FC1', 'C1': 'C2', 'C2': 'C1',
    'C3': 'C4', 'C4': 'C3', 'C5': 'C6', 'C6': 'C5', 'T7': 'T8', 'T8': 'T7', 'TP7': 'TP8', 'TP8': 'TP7',
    'CP5': 'CP6', 'CP6': 'CP5', 'CP3': 'CP4', 'CP4': 'CP3', 'CP1': 'CP2', 'CP2': 'CP1', 'P1': 'P2', 'P2': 'P1',
    'P3': 'P4', 'P4': 'P3', 'P5': 'P6', 'P6': 'P5', 'P7': 'P8', 'P8': 'P7', 'P9': 'P10', 'P10': 'P9',
    'PO7': 'PO8', 'PO8': 'PO7', 'PO3': 'PO4', 'PO4': 'PO3', 'O1': 'O2', 'O2': 'O1',
    'Fpz': 'Fpz', 'AFz': 'AFz', 'Fz': 'Fz', 'FCz': 'FCz', 'Cz': 'Cz', 'CPz': 'CPz', 'Pz': 'Pz', 'POz': 'POz', 'Oz': 'Oz'
}


# Create list of flipped
flipped_ch_names = []
flipped_indices = []
for ch in ch_names:
    flipped_ch = flipping_dict.get(ch, ch)
    flipped_ch_names.append(flipped_ch)
    flipped_indices.append(np.where(ch_names == flipping_dict.get(ch, ch))[0][0])

# DataFrame
global_electrode_idx = pd.DataFrame({
    'orig channel': ch_names,
    'original': np.arange(len(ch_names)),
    'flipped': flipped_indices,
    'flip channel': flipped_ch_names,
})

idx_original_electrodes = global_electrode_idx['original'].values
idx_flipped_electrodes = global_electrode_idx['flipped'].values

C3 = np.where((np.array(ch_names) == 'C3'))[0][0]
C4 = np.where((np.array((ch_names)) == 'C4'))[0][0]



# Long format with trial in line and variable in column (subject, trial, bin, stimulus, coherence, channel, bin)
# We divide analysis in 4 data represent 4 periodes (start at stimulus onset, stimulus offset, RS, and response), with following columns:

# | Subject | Trial | Bin | Stimulus | Coherence | Accuracy | TT | Response_moment | Value |
# | 1       | 1     | 1   | 0        | 0         | 1        | 5  | 1               | ...   |
# | 1       | 1     | 2   | 0        | 0         | 1        | 5  | 1               | ...   |
# | 1       | 1     | 3   | 0        | 0         | 1        | 5  | 1               | ...   |
# | 1       | 1     | ... | 0        | 0         | 1        | 5  | 1               | ...   |
# | 1       | 2     | 1   | 1        | 2         | 1        | ...| ..              | ...   |
# | 1       | 2     | 2   | 1        | 2         | ...      | .. | ...             | ...   |
# | 1       | 2     | 3   | 1        | 2         | ...      | .. | ...             | ...   |
# | ...     | ...   | ... | ...      | ...       | ...      | ...| ..              | ...   |
# | 22      | 960   | 10  | 1        | ...       | ...      | ...| ..              | ...   |


# Loop subjects, mean per conditions
subjects = np.arange(1,22)


x = 0
for s in tqdm(subjects, desc="Mean per conditions for subject", position=0, leave=False):

    # LOADING ----------------------------------------------------------------------
    #  Load epochs segmented for one subject ---
    epochs_onset  = np.load(path_eeg + f'\s{s}_stim_CSD_STFT.npy') # (960 trials, 64 eeg, X freq, X step time)
    epochs_offset = np.load(path_eeg + f'\s{s}_soff_CSD_STFT.npy')
    epochs_rs     = np.load(path_eeg + f'\s{s}_RS_CSD_STFT.npy')
    epochs_resp   = np.load(path_eeg + f'\s{s}_resp_CSD_STFT.npy')
    
    # Sanity check
    assert time_steps_onset.shape[0] == epochs_onset.shape[-1], "Time steps do not match with epochs SON data"
    assert time_steps_offset.shape[0] == epochs_offset.shape[-1], "Time steps do not match with epochs SOFF data"
    assert time_steps_rs.shape[0] == epochs_rs.shape[-1], "Time steps do not match with epochs RS data"
    assert time_steps_resp.shape[0] == epochs_resp.shape[-1], "Time steps do not match with epochs resp data"

    # Load behav data --- 
    data_behav = np.load(path_behav + f"\s{s}.npy") 
    trigg_resp = data_behav[:, 14].astype(np.float64).astype(int) # response trigger
    
    # From data export...
    data_suj = df_export[(df_export[:, 0] == s)] #extract data from subject
    RT = data_suj[:, 4]
    PMT = data_suj[:, 5]
    MT = data_suj[:, 6]
    stimLag = data_suj[:, 14]
    coherence = data_suj[:, 2]
    accuracy = data_suj[:, 3]
    response_moment = data_suj[:, 15]
    
    # Load EMG artifact vector data by subject -----
    vector_TT = np.load(path_tt + "\s" + str(s) + "_vector_TT_EMG_EEG.npy")
        # (0 Response before Response Signal (RS))
        # (1 Missed Trial)
        # (2 Anticipation RT < 200 ms) 
        # (3 Could not detect EMG onsets)
        # (4 technical problem) 
        # (5 Valid trial, to include)
        # (6 if EEG artifacts in valid trials)
    


    ######################### Reorganised electrod for MBL (all press is considere right = contro) #########################
    # find button press and recode variable 0 is left and 1 is right
    hand_press = []
    for t in range(len(trigg_resp)):
        if (trigg_resp[t] in [100, 101, 200,201]): # left
            hand_press.append(0)
        elif (trigg_resp[t] in [102, 103, 202,203]): # right
            hand_press.append(1)
        else:
            hand_press.append(-1)
    hand_press = np.array(hand_press)
    
    
    # === Flip electrode ipsi to contro for same hemisphere motor action ===
    
    # -- reorganize electrodes for contralateral response --
    EEG_stim = epochs_onset[:, :, :] # shape (960 trials, 64 eeg, frqs, x step epoch)
    EEG_soff = epochs_offset[:, :, :]
    EEG_rs = epochs_rs[:, :, :] 
    EEG_resp = epochs_resp[:, :, :] 

    EEG_stim_reorg = np.zeros_like(EEG_stim) 
    EEG_soff_reorg = np.zeros_like(EEG_soff)
    EEG_rs_reorg   = np.zeros_like(EEG_rs)
    EEG_resp_reorg = np.zeros_like(EEG_resp)
    
    for j in range (len(trigg_resp)):    
        
        # ---Create and add value eeg reorga if right response  ---
        if (trigg_resp[j] == 103):#correct right response, do not change anything
            EEG_stim_reorg[j,:,:] = EEG_stim[j,:,:]
            EEG_soff_reorg[j,:,:] = EEG_soff[j,:,:]
            EEG_rs_reorg[j,:,:]   = EEG_rs[j,:,:]
            EEG_resp_reorg[j,:,:] = EEG_resp[j,:,:]

        if (trigg_resp[j] == 101): #correct left response, ight hemisph activate ipsi, so flip for leavee all motor activation on same motor hemi with right hand = flip electrodes
            EEG_stim_reorg[j,:,:] = EEG_stim[j,idx_flipped_electrodes,:]
            EEG_soff_reorg[j,:,:] = EEG_soff[j,idx_flipped_electrodes,:]
            EEG_rs_reorg[j,:,:]   = EEG_rs[j,idx_flipped_electrodes,:]
            EEG_resp_reorg[j,:,:] = EEG_resp[j,idx_flipped_electrodes,:] 
 
        if (trigg_resp[j] == 102):#incorrect right response. Expected is left. 
            EEG_stim_reorg[j,:,:] = EEG_stim[j,:,:]
            EEG_soff_reorg[j,:,:] = EEG_soff[j,:,:]
            EEG_rs_reorg[j,:,:]   = EEG_rs[j,:,:]
            EEG_resp_reorg[j,:,:] = EEG_resp[j,:,:]   

        if (trigg_resp[j] == 100):#incorrect left response, right hemisph activate ipsi, so flip for leavee all motor activation on same motor hemi with right hand
            EEG_stim_reorg[j,:,:] = EEG_stim[j,idx_flipped_electrodes,:]
            EEG_soff_reorg[j,:,:] = EEG_soff[j,idx_flipped_electrodes,:]
            EEG_rs_reorg[j,:,:]   = EEG_rs[j,idx_flipped_electrodes,:]
            EEG_resp_reorg[j,:,:] = EEG_resp[j,idx_flipped_electrodes,:]


    
    # ------- Select channels of interest for MBL -------
    # Note.   C3 = 12 (is left electrod so right pres). In flip data, C3 is controlateral rightto reponse
            # C4 = 49 (is right electrod so left press, In flip data, C4 is ipsilatral to right response
            # As twomey et al (2016) : "Based on inspection of this topography, we chose the standard left- and right-hemisphere motor sites C3 and C4 
            # for -band amplitude measurements. Stimulus- and response-aligned -band lateralization was calculated by subtracting amplitude values 
            # at the electrode ipsilateral to the responding hand from those at the electrode contralateral to the responding hand (contralateral ipsilateral)."
    
    # Select MB freq 8-30Hz
    frqs = np.where(freq_stft <= 30)[0]
    
    MBL_stim = np.mean(EEG_stim_reorg[:, C3, frqs, :] - EEG_stim_reorg[:, C4, frqs, :], axis = 1)   # in mV/m² -- (960 trials, x frqs, x step epoch)
    MBL_soff = np.mean(EEG_soff_reorg[:, C3, frqs, :] - EEG_soff_reorg[:, C4, frqs, :] , axis = 1)
    MBL_rs = np.mean(EEG_rs_reorg[:, C3, frqs,:] - EEG_rs_reorg[:, C4, frqs,:], axis = 1)
    MBL_resp = np.mean(EEG_resp_reorg[:, C3, frqs,:] - EEG_resp_reorg[:, C4, frqs,:] , axis = 1)


    # Ordonnate trials in columns for long data -------
    # -> columns: Subject | Trial | Bin | Stimulus | Coherence | Accuracy | TT | Response_moment | Value
    
    for trial in range(MBL_resp.shape[0]): # 960 trials max
        # Find value voltage of corresponding bin --
        # Remember epoch are aligned (and centered at 0) on corresponding event (stimulus onset, stimulus offset, RS, resp)
        
        # ------ Stimulus period ------
        # Select coresponding bin (-0 correspond time_TF_onset[15])
        idx_start_bin = index_nearest(time_steps_onset, bin_stim[0]) # index of stimulus onset in epoch time
        idx_end_bin = index_nearest(time_steps_onset, bin_stim[-1])
        epoch_trial = MBL_stim[trial, idx_start_bin:idx_end_bin] # from -0.1 to 2s (stimulus duration max)
        

        for b in range(len(bin_stim) -1): # -1 because we take bin to bin+1
            temp_voltage_bin = epoch_trial[b]

            # Get and add trial info to data container       
            arr_stim[row_idx_stim, :] = [
                x+1,
                trial+1,
                b, # bin 0 is -0.1s to 0s average bin
                stimLag[trial],
                coherence[trial],
                accuracy[trial],
                vector_TT[trial],
                response_moment[trial],
                temp_voltage_bin
            ]
            row_idx_stim += 1 # next row in data container (same trial, next bin)


        # ------ Delayed period ------
        idx_start_bin = index_nearest(time_steps_offset, bin_delay[0]) 
        idx_end_bin = index_nearest(time_steps_offset, bin_delay[-1])
        epoch_trial = MBL_soff[trial, idx_start_bin:idx_end_bin]
        

        for b in range(len(bin_delay) -1):
            temp_voltage_bin = epoch_trial[b]
        
            # Get and add trial info to data container       
            arr_delay[row_idx_delay, :] = [
                x+1,
                trial+1,
                b,
                stimLag[trial],
                coherence[trial],
                accuracy[trial],
                vector_TT[trial],
                response_moment[trial],
                temp_voltage_bin
            ]
            row_idx_delay += 1
    
    
        #  RS period ------
        idx_start_bin = index_nearest(time_steps_rs, bin_rs[0]) 
        idx_end_bin = index_nearest(time_steps_rs, bin_rs[-1])
        epoch_trial = MBL_rs[trial, idx_start_bin:idx_end_bin]

        for b in range(len(bin_rs) -1):
            temp_voltage_bin = epoch_trial[b]
        
            # Get and add trial info to data container       
            arr_rs[row_idx_rs, :] = [
                x+1,
                trial+1,
                b,
                stimLag[trial],
                coherence[trial],
                accuracy[trial],
                vector_TT[trial],
                response_moment[trial],
                temp_voltage_bin
            ]
            row_idx_rs += 1 # next row in data container (same trial, next bin)
    
    
        # Response period ------
        idx_start_bin = index_nearest(time_steps_resp, bin_resp[0]) 
        idx_end_bin = index_nearest(time_steps_resp, bin_resp[-1])
        epoch_trial = MBL_resp[trial, idx_start_bin:idx_end_bin]

        for b in range(len(bin_resp) -1):
            temp_voltage_bin = epoch_trial[b]
        
            # Get and add trial info to data container       
            arr_resp[row_idx_resp, :] = [
                x+1,
                trial+1,
                b,
                stimLag[trial],
                coherence[trial],
                accuracy[trial],
                vector_TT[trial],
                response_moment[trial],
                temp_voltage_bin
            ]
            row_idx_resp += 1 # next row in data container (same trial, next bin)
    
    
    # Average LRP per conditions also ---------------------------------------------------------------------
    for stim_duration in range(2): 
        for coh in range(3):
            for press in range(3):
                if (press == 2): # all button pressed is controlateral to response on same hemisphere (right press, because we flip electrods)
                    mask = (vector_TT == 5)  & (stimLag == possible_stim[stim_duration]) & (coherence == possible_coh[coh]) & (response_moment == 1) & (accuracy != -1)
                else: # differs button press (right and lef, to plot topo as Kelly et al. 2013)
                    mask = (vector_TT == 5)  & (stimLag == possible_stim[stim_duration]) & (coherence == possible_coh[coh]) & (response_moment == 1)  & (accuracy != -1) & (hand_press == press)
            
                # == For all period, mean trials per conditions ==
                #  Stim ON 
                MB_onset_conditions[x, stim_duration, coh,:,:, :, press] = np.nanmean(EEG_stim_reorg[mask], axis = 0)
                MB_onset_conditions_no_reorg[x, stim_duration, coh,:,:, :,press] = np.nanmean(EEG_stim[mask], axis = 0) 
                
                # Stim OFF
                MB_offset_conditions[x, stim_duration, coh,:,:, :, press] = np.nanmean(EEG_soff_reorg[mask], axis = 0)
                MB_offset_conditions_no_reorg[x, stim_duration, coh,:,:, :, press] = np.nanmean(EEG_soff[mask], axis = 0)
                
                #  RS 
                MB_rs_conditions[x, stim_duration, coh,:,:, :, press] = np.nanmean(EEG_rs_reorg[mask], axis = 0)
                MB_rs_conditions_no_reorg[x, stim_duration, coh,:,:, :, press] = np.nanmean(EEG_rs[mask], axis = 0) 

                #  resp 
                MB_resp_conditions[x, stim_duration, coh,:,:,:, press] = np.nanmean(EEG_resp_reorg[mask], axis = 0)
                MB_resp_conditions_no_reorg[x, stim_duration, coh,:,:,:, press] = np.nanmean(EEG_resp[mask], axis = 0) 
    
    x += 1
print("--- End of loop for all subjects ---")


os.chdir(current)
# # Save average data by condition for plot -----------------------------
# # --------- reorga
# np.save(f"av_epochs_onset_reorg.npy", MB_onset_conditions) # (all subjects, 2 stim, 3coh, 64 eeg, x freqs, time_steps)
# np.save(f"av_epochs_offset_reorg.npy", MB_offset_conditions)
# np.save(f"av_epochs_rs_reorg.npy", MB_rs_conditions) 
# np.save(f"av_epochs_resp_reorg.npy", MB_resp_conditions) 

# # -------- no reorga 
# np.save(f"av_epochs_onset_no_reorg.npy", MB_onset_conditions_no_reorg)
# np.save(f"av_epochs_offset_no_reorg.npy", MB_offset_conditions_no_reorg)
# np.save(f"av_epochs_rs_no_reorg.npy", MB_rs_conditions_no_reorg)
# np.save(f"av_epochs_resp_no_reorg.npy", MB_resp_conditions_no_reorg)


# Save long data format ... ------------------------------------------------
# Delete unused rows
df_arr_stim = arr_stim[arr_stim[:, 0] != 0 ]
df_arr_delay = arr_delay[arr_delay[:, 0] != 0 ]
df_arr_rs = arr_rs[arr_rs[:, 0] != 0 ]
df_arr_resp = arr_resp[arr_resp[:, 0] != 0 ]
    
# long format DataFrame
df_long_stim = pd.DataFrame(df_arr_stim, columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'accuracy', 'TT', 'response_moment', 'mean_voltage'])
df_long_delay = pd.DataFrame(df_arr_delay, columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'accuracy', 'TT', 'response_moment', 'mean_voltage'])
df_long_rs= pd.DataFrame(df_arr_rs, columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'accuracy', 'TT', 'response_moment', 'mean_voltage'])
df_long_resp= pd.DataFrame(df_arr_resp, columns=['subject', 'trial', 'bin', 'stimulus', 'coherence',  'accuracy', 'TT', 'response_moment', 'mean_voltage'])
df_long_resp.loc[(df_long_resp["TT"] == 5) & (df_long_resp["mean_voltage"].isna()), "TT"] = 4  # convert 16 NA in TT 4 (if abnormal pretreatment) 


# # Save PMA data
# os.chdir(current)
# df_long_stim.to_csv('df_long_stim_MBL.csv', index=False) 
# df_long_delay.to_csv('df_long_delay_MBL.csv', index=False) 
# df_long_rs.to_csv('df_long_rs_MBL.csv', index=False) 
# df_long_resp.to_csv('df_long_resp_MBL.csv', index=False) 
    
# df_long_stim  = pd.read_csv('df_long_stim_MBL.csv')
# df_long_delay = pd.read_csv('df_long_delay_MBL.csv')
# df_long_rs    = pd.read_csv('df_long_rs_MBL.csv')
# df_long_resp  = pd.read_csv('df_long_resp_MBL.csv')
    
    
    
# ...and convert long format to wide format for JASP analysis ---------------
def prepare_jasp_csv(df, bin_min=0, bin_max= None, accuracy_exclude=-1):
    
    #  Filtrage ------------------
    df_filtered = df[
        (df['TT'] == 5) &
        (df['response_moment'] == 1) &
        (df['accuracy'] != accuracy_exclude) &
        (df['bin'] >= bin_min) &
        (df['bin'] <= bin_max)
    ].copy()
    
    #  Convertir en variables catégorielles ------------------
    df_filtered['subject'] = df_filtered['subject'].astype('category')
    df_filtered['trial'] = df_filtered['trial'].astype('category')

    #  Centrer les variables ------------------
    df_filtered['coherence_norm'] = df_filtered['coherence'] - df_filtered['coherence'].mean()
    df_filtered['stimulus_norm'] = df_filtered['stimulus'] - df_filtered['stimulus'].mean()
    df_filtered['mean_voltage'] = df_filtered['mean_voltage']/1000  # si pas de centrage requis - V/m² to mV/m²

    # # Vérification des NAs et plage ------------------
    # print(df_filtered.isna().sum())
    # print(df_filtered['mean_voltage'].min(), df_filtered['mean_voltage'].max())


    #  Calculer les moyennes par condition pour RM-ANOVA ------------------
    df_means = df_filtered.groupby(
        ['subject', 'stimulus', 'coherence', 'bin'],
        as_index=False
    )['mean_voltage'].mean()

    #  Format wide pour JASP ------------------
    df_means['condition'] = (
        'stim' + df_means['stimulus'].astype(int).astype(str) +
        '_coh' + (df_means['coherence']*100).astype(int).astype(str) +
        '_bin' + df_means['bin'].astype(int).astype(str)
    )

    #  Moyenne par sujet × condition pour pivot ------------------
    df_means_subject = df_means.groupby(
        ['subject', 'condition'],
        as_index=False
    )['mean_voltage'].mean()
    
    df_wide = df_means_subject.pivot(index='subject', columns='condition', values='mean_voltage').reset_index()
    
    
    def sort_key(col_name):
        match = re.match(r'stim(\d+)_coh(\d+)_bin(\d+)', col_name)
        if match:
            return tuple(int(x) for x in match.groups())
        else:
            return (0,0,0)  # pour 'subject' ou colonnes inattendues

    cols = df_wide.columns.tolist()
    cols_sorted = [cols[0]] + sorted(cols[1:], key=sort_key)  # conserver 'subject' en premier
    df_wide = df_wide[cols_sorted]
    
    return df_wide


# # Stimulus period ------------------
# df_wide = prepare_jasp_csv(df_long_stim, bin_min=0, bin_max=20)
# df_wide.to_csv("MBL_stim_JASP.csv.csv", index=False)

# #  Delay period ------------------
# df_wide = prepare_jasp_csv(df_long_delay, bin_min=10, bin_max=30) # 1sec after stim offset
# df_wide.to_csv("MBL_delay_JASP.csv", index=False)

# # RS period ------------------
# df_wide = prepare_jasp_csv(df_long_rs, bin_min=1, bin_max= 3) # 1st after RS onset (+300ms)
# df_wide.to_csv("MBL_rs_JASP.csv", index=False)

# # Response locked period ------------------
# df_wide = prepare_jasp_csv(df_long_resp, bin_min=1, bin_max= 3) 
# df_wide.to_csv("MBL_resp_JASP_bis.csv", index=False)


print('MBL average by condition processed.')
# [X]