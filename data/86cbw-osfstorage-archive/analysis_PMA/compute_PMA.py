""" 
Find number of PMA per 0.1 bin with EMG segmented 

---
Note. Partial muscle activations (PMAs) are subthreshold EMG bursts that do not
produce an overt response. According to the extended evidence accumulation
theory, PMAs occur when noisy motor preparatory activity transiently exceeds
gating inhibition. To estimate their time-varying occurrence rate, PMAs were
averaged within 0.1 s time bins separately for each subject, EMG channel,
and experimental condition.
"""

import numpy as np
import pandas as pd
import os
import copy
import sys
import re
import matplotlib.pyplot as plt
current = os.path.join(os.getcwd()) 
toolbox_EMG_folder = r"C:\Users\erc\iCloudDrive\ACHARD_these\EMG_method\toolbox_EMG"
sys.path.extend([toolbox_EMG_folder]) # to add toolbox path to environment paths
import emgTools as emg
import evtTools as evt
from Events import Events
import mne

datapath_behav = os.path.abspath(os.path.join(current, "..", "..", "data", "behavioral"))
datapath_EMG_segmented = os.path.abspath(os.path.join(current, "..", "..", "data", "emg_segmented"))
datapath_vector_TT = os.path.abspath(os.path.join(current, "..", "..", "data", "vector_trial_type"))
datapath_JASP = current # jasp saving


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


# ------------------------------------------------ ID triggers ------------------------------------------------
"""
Define triggers index of experiment
"""

emg_channels_idx = {0: 'EMGL', 1: 'EMGR'}
trigg_stim =  [21, 22, 23, 24, 25, 26,     31, 32, 33, 34, 35, 36]
trigg_stim_id = {'stim_left_300_0':21,'stim_left_300_12':22,'stim_left_300_50':23,
                'stim_left_2000_0':31,'stim_left_2000_12':32,'stim_left_2000_50':33,

                'stim_right_300_0':24,'stim_right_300_12':25,'stim_right_300_50':26,
                'stim_right_2000_0':34,'stim_right_2000_12':35,'stim_right_2000_50':36} # used for segmentation    

trigg_RS =   [121, 122, 123, 124, 125, 126,     131, 132, 133, 134, 135, 136]  
trigg_RS_id = {'stim_left_300_2':121,'stim_left_300_11':122,'stim_left_300_40':123,
                'stim_left_2000_2':131,'stim_left_2000_11':132,'stim_left_2000_40':133,

                'stim_right_300_2':124,'stim_right_300_11':125,'stim_right_300_40':126,
                'stim_right_2000_2':134,'stim_right_2000_11':135,'stim_right_2000_40':136}  # used for segmentation for response signal (trigg stimulys +100)

trigg_resp_unique = 101,102,103,104
trigg_resp_id = {'left_correct': 101, 'left_incorrect': 100, 'right_correct' : 103, 'right_incorrect':102, 
                    'left_correct_beforeRS': 201, 'left_incorrect_beforeRS': 200, 'right_correct_beforeRS' : 203, 'right_incorrect_beforeRS':202, 
                    'no_resp':250}

trigg_emg_id = {'onset': 350 , 'offset': 351 }





# ------------------------------------------------ Data containers ------------------------------------------------
"""
Create data containers for PMA
- for all trials
- per period (stimulus, delay, response)
"""

# Bin length (in sample) here 0.1s --- 
# e.g., 50ms is 0.05s * 1024Hz = 51.2 samples
# e.g., 100ms is 0.1s * 1024Hz = 102.4 samples
len_bin = 0.1*1024  

# Create bin ---
max_after_RS = 0.3
bin_all = np.arange(len_bin, (2+3+ max_after_RS) * 1024 , len_bin)  # from stimulus onset 0 to +0.2s after RS onset
bin_stim = np.arange(len_bin, (2) * 1024 +1, len_bin)      # max 2s stimulus duration
bin_delay = np.arange(len_bin, (3) * 1024 +1 , len_bin)     # 3s delay
bin_resp = np.arange(len_bin, (max_after_RS) * 1024 +1, len_bin)  # 200ms after RS

# Data containers for PMA in each bin ---
data_bin_all = np.zeros((22, 960, 2,3,2, len(bin_all)))  # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)

data_bin_stim = np.zeros((22, 960,  2,3,2, len(bin_stim)))  # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
data_bin_delay = np.zeros((22, 960, 2,3,2, len(bin_delay)))  
data_bin_resp = np.zeros((22, 960,  2,3,2, len(bin_resp)))  
data_bin_resp_lock = np.zeros((22, 960,  2,3,2, len(bin_resp))) 

# Check sample and timing of experiment
data_len_sample_delay = np.full((22, 960), np.nan)
pma_after_last_bin = np.zeros((22,))  # count number of PMA after last bin in delay period (due to frame lag)
len_off_rs = np.zeros((22, 960))  # length of delay period in sample (from offset to RS onset, expected 3072 samples = 3s)


# Data containers for long format ---
n_total = 22 * 960 * len(bin_stim) * 2   # subject * trials * bin * channels
arr_stim = np.zeros((n_total, 10))       # 10 cols: Subject, Trial, Bin, Stimulus, Coherence, Channel, Accuracy, TT, Response_moment, PMA
row_idx = 0

n_total = 22 * 960 * len(bin_delay) * 2 
arr_delay = np.zeros((n_total, 10))
row_idx_delay = 0

n_total = 22 * 960 * len(bin_resp) * 2 
arr_resp = np.zeros((n_total, 10))
row_idx_resp = 0

n_total = 22 * 960 * len(bin_resp) * 2 
arr_resp_onset =  np.zeros((n_total, 10))
row_idx_resp_onset = 0


# Possible level ---
possible_coh_levels =  [0, 0.12, 0.5]
possible_stim_levels = [18, 120] # 0.3s / 2s
possible_channel =     [0, 1]    # opp / same 
    



# ------------------------------------------------ Compute PMA ------------------------------------------------

s = 0
for x in range(22):
    nameSubj = 's' + str(x+1)

    # ------------------------ Load data -----------------------------------
    # Load individual data
    os.chdir(datapath_behav)
    data_behav = np.load(nameSubj +".npy") 

    coherence = data_behav[:, 5].astype(np.float64)
    stim_pres = data_behav[:, 6].astype(np.float64).astype(np.int64) 
    response_moment = data_behav[:,10].astype(np.float64).astype(np.int64)
    acc = data_behav[:, 11].astype(np.float64).astype(int) 
    rt = data_behav[:, 12].astype(np.float64)
    motion_direction = data_behav[:, 7].astype(np.float64).astype(np.int64)  # 180: left   0: right 
    key_pressed = data_behav[:,9].astype(np.float64).astype(np.int64) 
    trigg_resp = data_behav[:, 14].astype(np.float64).astype(int)


    # Load EMG artifact vector data by subject
    vector_TT = np.load(datapath_vector_TT + "\s" + str(x+1) + "_vector_TT_EMG_EEG.npy")
        # (0 Response before Response Signal)
        # (1 Missed Trial)
        # (2 Anticipation RT < 200 ms) 
        # (3 Could not detect EMG onsets)
        # (4 technical problem) 
        # (5 Valid trial, to include)
        # (6 if EEG artifacts)

        
    # Load segmented events
    sf = 1024.0
    folder = datapath_EMG_segmented
    MARK_seg = evt.loadSegmented(os.path.join(folder,nameSubj + "_segmentedEvts.csv"),headerlines=1, sep=',', sf=sf,
                        colSample=0, colTime=1, colCode=2, colChan=3,
                        colTrialIdx=4, 
                        colTminSample=5, colTminTime=6,
                        colT0Sample=7, colT0Time=8, 
                        colTmaxSample=9, colTmaxTime=10)
    MARK_seg_modified = copy.deepcopy(MARK_seg)

    

    # -------------------------- Trial {i} ----------------------------
    for i in range(MARK_seg.nbTrials()):
        
        # Exclude trials with Biosemi technical problem (add 99 if few trial)
        if (MARK_seg.listEvtsTrials[i].code[-2] == 99): # (for auto segmented) or (s == 3 and i ==106 or i == 435 or i == 562) or (s==7 and i== 806) or (s==9 and i== 616) or (s==14 and i== 564 or i ==713 or i ==714)  or (s==16  and i == 350 or i ==936) or (s==19  and i == 699) or (s==20  and i == 273 or i ==741 or i ==814):  
            data_bin_all[s, i,:, :, :, :] = np.nan
            len_off_rs[s, i] = np.nan
            
        # Exclude no valid trials (TT != 5) 
        elif (vector_TT[i] != 5): 
            data_bin_all[s, i,:, :, :, :] = np.nan
            len_off_rs[s, i] = np.nan
        
        # Find PMA in valid trial (and add check for accuracy (correct+incorrect) and response moment (== 1) is valid)
        elif (vector_TT[i] == 5) and (acc[i] != -1) and (response_moment[i] == 1):
            
            # == 1 == Find channel_associated_with_response is the channel of the response --------------------------------------
            if trigg_resp[i] in [100,101]:
                channel_associated_with_response = 0 # for incorrect/correct left response channel of the response is left (0)
                channel_not_associated_with_response = 1
            elif trigg_resp[i] in [102,103]:
                channel_associated_with_response = 1  # for incorrect/correct right response channel of the response is right (1)
                channel_not_associated_with_response = 0    



            # == 2 == Find stim ON, OFF, RS, PMT burst and resp trigger  -------------------------------------------------------------------
            # Here, we find index of interest triggers: stimulus offset, response signal, EMG onset of response, and response
            idx_onset  = int(MARK_seg.listEvtsTrials[i].findEvents(code=list(trigg_stim_id.values()))) 
            idx_offset = int(MARK_seg.listEvtsTrials[i].findEvents(code=1))
            idx_RS     = int(MARK_seg.listEvtsTrials[i].findEvents(code=list(trigg_RS_id.values()))) 
            idx_resp   = int(MARK_seg.listEvtsTrials[i].findEvents(code=trigg_resp[i])) 

            # for EMG onset of response, we start from response to RS to find the first 350 (correspond to onset)
            nb_trigg= np.flip(np.arange(idx_resp)) # reverse trigger to check EMG onset before response. So first trigger is always 350 (the two next are response and 351)
            for n in nb_trigg:
                if  MARK_seg.listEvtsTrials[i].code[n] == 350:
                    if (MARK_seg.listEvtsTrials[i].chan[n] == channel_associated_with_response): 
                        check = 1
                        idx_pmt_burst = n # save position of last burst of PMT
                        break



            # == 3 == Create bin of sample, check if PMA inside (+1) -------------------------------------------------------------------------
            # Find variable of interest of trial
            stim_pres_trial = stim_pres[i]
            coh_trial = coherence[i]
            rt_trial = rt[i]
            sample_trial = MARK_seg.listEvtsTrials[i].Latency.sample     
            code_trial = MARK_seg.listEvtsTrials[i].code
            channel_trial = MARK_seg.listEvtsTrials[i].chan
            latence_trial = MARK_seg.listEvtsTrials[i].Latency.time  
            
            # -- Note--. Be careful stimulus onset start at +300ms after start of trial (correspond to baseline -0.3s to 0)
            # In consequence, we realign 0 locked on stimulus onset to create bin with idx_onset of trials
            sample_trial = sample_trial - sample_trial[idx_onset]
            
            
            # Index of conditions (stimulus and coherence)
            if  stim_pres_trial == 18:  # 300 ms
                stim = 0 
                no_stim = 1
            elif stim_pres_trial == 120: # 2s
                stim = 1 
                no_stim = 0  
            
            if   coh_trial == 0:    coh = 0
            elif coh_trial == 0.12: coh = 1
            elif coh_trial == 0.5:  coh = 2  
            

            
            # == 3a == From stimulus onset to 0.2s after RS --------------------------------------------------------------
            for idx in range(idx_onset, idx_pmt_burst): 
                # Find informations of correspond trigger in event
                temp_sample = sample_trial[idx]
                temp_code   = code_trial[idx]
                temp_chan   = channel_trial[idx]
                temp_lat    = latence_trial[idx] 
                
                # -- Navigate in bin (previous create) 
                for b in range(len(bin_all)):
                    # Check if actual sample is in bin ...
                    if (bin_all[b] - len_bin) < temp_sample <= bin_all[b]:  

                        # ... and if corresponding code is an EMG onset (350) ..                  
                        if temp_code == 350:
                            # .. if yes add +1 PMA to corresponding bin and channel (same or opposite as response)
                            if temp_chan == channel_associated_with_response:
                                data_bin_all[s, i, stim, coh, 1, b] += 1    # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
                                
                            elif temp_chan == channel_not_associated_with_response:
                                data_bin_all[s, i, stim, coh, 0, b] += 1
                    
                    # find idx_pmt_burst, and add nan for rest of bin (so after EMG onset of response include)
                    elif sample_trial[idx_pmt_burst] <= (bin_all[b] - len_bin):
                        # # Debug...
                        # print(f'{sample_trial[idx_pmt_burst]}  <= {(bin_all[b] - len_bin)}')
                        data_bin_all[s, i, stim, coh, :, b:] = np.nan
                        break

            
            
            
            # == 3b == From stimulus, delay and response period locked -----------------------------------------------------
            
            # ---- (A) Stimulus period (stimulus onset locked) -------------------------------------------------------------
            sample_trial = MARK_seg.listEvtsTrials[i].Latency.sample   
            # Note. Here we realign 0 to stimulus onset (start of stimulus period)
            sample_trial = sample_trial - sample_trial[idx_onset]
            
            # for plot to see after offstet for 300ms
            if stim_pres_trial == 120 : 
                idx_stim_plot = idx_offset  
            elif stim_pres_trial == 18:
                idx_stim_plot = idx_RS  # +3s after offset but break wheen 
            
            for idx in range(idx_onset, idx_stim_plot):
                # Find informations of correspond trigger in event
                temp_sample = sample_trial[idx]
                temp_code   = code_trial[idx]
                temp_chan   = channel_trial[idx]
                temp_lat    = latence_trial[idx] 
                
                # Stop after 2s (to plot acccumation after for 300ms condition)
                if temp_sample > (2*1024) : # in Hz 
                    break
                
                # -- Navigate in bin (previous create) 
                for b in range(len(bin_stim)):
                    # Check if actual sample is in bin ...
                    if (bin_stim[b] - len_bin) < temp_sample <= bin_stim[b]:    
                        # ... and if corresponding code is an EMG onset (350) ..                  
                        if temp_code == 350:
                            # .. if yes add +1 PMA to corresponding bin and channel (same or opposite as response)
                            if temp_chan == channel_associated_with_response:
                                data_bin_stim[s, i, stim, coh, 1, b] += 1    # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
                            elif temp_chan == channel_not_associated_with_response:
                                data_bin_stim[s, i, stim, coh, 0, b] += 1
                    
                        # PMA after last theorical sample (different due to frame lag)
                        if temp_sample > bin_stim[b]:
                            print("PMA no in last bin [stim]")
            
            # Add line for long format
            for b in range(len(bin_stim)):
                for ch in range(2):
                    arr_stim[row_idx, :] = [s+1, i+1, b, stim, coh, ch, acc[i], vector_TT[i], response_moment[i], data_bin_stim[s, i, stim, coh, ch, b]]
                    row_idx += 1
            
            
            
            # ---- (B) Delay period (stimulus offset locked) ----------------------------------------------------------------
            sample_trial = MARK_seg.listEvtsTrials[i].Latency.sample  
            sample_trial = sample_trial - sample_trial[idx_offset]    # Note. Here we realign 0 to stimulus offset (start of delay period)
            
            # -- Compute length of bin in dealy period 
            len_off_rs[s, i] = sample_trial[idx_RS] - sample_trial[idx_offset]
            
            for idx in range(idx_offset, idx_RS):
                # Find informations of correspond trigger in event
                temp_sample = sample_trial[idx]
                temp_code   = code_trial[idx]
                temp_chan   = channel_trial[idx]
                temp_lat    = latence_trial[idx] 
                
                # -- Navigate in bin (previous create) 
                # Note. Last bins can be varied due to PsychoPy frame lag. So RS apear exactly at the 180e frame after offset. For few trials +/- 1 or 2 frames (0.016s) can be add.
                # in bin segmentation, bins are equal and the last also (100ms = 102.4 sample). If a frame is introduce after 3072 sample is compute here.
                
                for b in range(len(bin_delay)):
                    # Check if actual sample is in bin ...
                    if (bin_delay[b] - len_bin) < temp_sample <= bin_delay[b]:    
                        # ... and if corresponding code is an EMG onset (350) ..                  
                        if temp_code == 350:
                            # .. if yes add +1 PMA to corresponding bin and channel (same or opposite as response)
                            if temp_chan == channel_associated_with_response:
                                data_bin_delay[s, i, stim, coh, 1, b] += 1    # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
                            elif temp_chan == channel_not_associated_with_response:
                                data_bin_delay[s, i, stim, coh, 0, b] += 1
                            
                        # PMA after last theorical sample 3072 correspond to 3s delay (different due to frame lag)
                        if temp_sample > bin_delay[b]:
                            print("PMA no in last bin [delay]")
                            pma_after_last_bin[s] += 1
            
            # Add line for long format
            for b in range(len(bin_delay)):
                for ch in range(2):
                    arr_delay[row_idx_delay, :] = [s+1, i+1, b, stim, coh, ch, acc[i], vector_TT[i], response_moment[i], data_bin_delay[s, i, stim, coh, ch, b]]
                    row_idx_delay += 1
                                       
                                       
            
            # ---- (C) Response period (rs locked) -------------------------------------------------------------------------
            sample_trial = MARK_seg.listEvtsTrials[i].Latency.sample   
            # Note. Here we realign 0 to rs (start of response period)
            sample_trial = sample_trial - sample_trial[idx_RS]
            
            # We exclude EMG onset burst of response bevause it is not a PMA
            for idx in range(idx_RS, idx_pmt_burst):
                # Find informations of correspond trigger in event
                temp_sample = sample_trial[idx]
                temp_code   = code_trial[idx]
                temp_chan   = channel_trial[idx]
                temp_lat    = latence_trial[idx] 
                
                # -- Navigate in bin (previous create) 
                for b in range(len(bin_resp)):
                    # Check if actual sample is in bin ...
                    if (bin_resp[b] - len_bin) < temp_sample <= bin_resp[b]:    
                        # check if temp_sample is < idx_pmt_burst, else add np.nan to not take account bin after RT in response period
                        if (idx >= idx_pmt_burst-1):  # last loop add nan to rest of bins
                            data_bin_resp[s, i, stim, coh, :, b+1:] = np.nan
                            
                        # ... and if corresponding code is an EMG onset (350) ..                  
                        if temp_code == 350:
                            # .. if yes add +1 PMA to corresponding bin and channel (same or opposite as response)
                            if temp_chan == channel_associated_with_response:
                                data_bin_resp[s, i, stim, coh, 1, b] += 1    # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
                            elif temp_chan == channel_not_associated_with_response:
                                data_bin_resp[s, i, stim, coh, 0, b] += 1
            
            # Add line for long format
            for b in range(len(bin_resp)):
                for ch in range(2):
                    arr_resp[row_idx_resp, :] = [s+1, i+1, b, stim, coh, ch, acc[i], vector_TT[i], response_moment[i], data_bin_resp[s, i, stim, coh, ch, b]]
                    row_idx_resp += 1
                    
            
            # ---- (D) Response period (resp locked) -------------------------------------------------------------------------
            sample_trial = MARK_seg.listEvtsTrials[i].Latency.sample   
            # Note. Here we realign 0 on response (start of response period)
            sample_trial = sample_trial - sample_trial[idx_pmt_burst] #idx_pmt_burst
            
            nb_trigg = np.flip(np.arange(idx_RS, idx_pmt_burst)) # to PMT burst at RS
            bin_resp_lock = -bin_resp #-np.flip(bin_resp)
            
            # We exclude EMG onset burst of response because it is not a PMA (is EMG onset of response)
            for idx in nb_trigg:
                # Find informations of correspond trigger in event
                temp_sample = sample_trial[idx]
                temp_code   = code_trial[idx]
                temp_chan   = channel_trial[idx]
                temp_lat    = latence_trial[idx] 
                               
                # -- Navigate in bin (previous create) 
                for b in range(len(bin_resp_lock)):
                    # Check if actual sample is in bin ...
                    if (bin_resp_lock[b] - len_bin) < temp_sample <= bin_resp_lock[b]:    
                        
                        # check if temp_sample is < idx_pmt_burst, else add np.nan to not take account bin after RT in response period
                        if (idx >= idx_pmt_burst-1):  # last loop add nan to rest of bins
                            data_bin_resp_lock[s, i, stim, coh, :, b+1:] = 0
                            
                        # ... and if corresponding code is an EMG onset (350) ..                  
                        if temp_code == 350:
                            # .. if yes add +1 PMA to corresponding bin and channel (same or opposite as response)
                            if temp_chan == channel_associated_with_response:
                                data_bin_resp_lock[s, i, stim, coh, 1, b] += 1    # 22 subject, 960 trials, 2 stim, 3 coh, 2 chan, X bin)
                            elif temp_chan == channel_not_associated_with_response:
                                data_bin_resp_lock[s, i, stim, coh, 0, b] += 1

            # Add line for long format
            for b, bin_edge in enumerate(bin_resp_lock):
                for ch in range(2):
                    arr_resp_onset[row_idx_resp_onset, :] = [s+1, i+1, b, stim, coh, ch, acc[i], vector_TT[i], response_moment[i], data_bin_resp_lock[s, i, stim, coh, ch, b]]
                    row_idx_resp_onset += 1
    
    
    s+=1        
    print(f"Subject {x+1} completed")
print(f"    >>> Processing completed on {s} subjects")     


# Long format data ----
df_long_stim = pd.DataFrame(arr_stim[:row_idx, :], columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'channel', 'accuracy', 'TT', 'response_moment', 'PMA'])
df_long_delay = pd.DataFrame(arr_delay[:row_idx_delay, :], columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'channel', 'accuracy', 'TT', 'response_moment', 'PMA'])
df_long_resp = pd.DataFrame(arr_resp[:row_idx_resp, :], columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'channel', 'accuracy', 'TT', 'response_moment', 'PMA'])
df_long_resp_onset = pd.DataFrame(arr_resp_onset[:row_idx_resp, :], columns=['subject', 'trial', 'bin', 'stimulus', 'coherence', 'channel', 'accuracy', 'TT', 'response_moment', 'PMA'])

# replace per coherence real levels
df_long_stim['coherence'] = df_long_stim['coherence'].replace({0:0 ,1:0.12, 2:0.50})
df_long_delay['coherence'] = df_long_delay['coherence'].replace({0:0 ,1:0.12, 2:0.50})
df_long_resp['coherence'] = df_long_resp['coherence'].replace({0:0 ,1:0.12, 2:0.50})
df_long_resp_onset['coherence'] = df_long_resp_onset['coherence'].replace({0:0 ,1:0.12, 2:0.50})

# # Save PMA data long format if need
# df_long_stim.to_csv('df_long_stim_PMA.csv', index=False) 
# df_long_delay.to_csv('df_long_delay_PMA.csv', index=False) 
# df_long_resp.to_csv('df_long_resp_PMA.csv', index=False) 
# df_long_resp_onset.to_csv('df_long_resp_onset_PMA.csv', index=False) 








# =============================================================== Average trials by conditions ===============================================================
# Create data containers
av_data_bin_all   = np.zeros((22, 2,3,2, len(bin_all)))  # 22 subjects,  2 stim, 3 coh, 2 chan, X bin)
av_data_bin_stim  = np.zeros((22, 2,3,2, len(bin_stim))) 
av_data_bin_delay = np.zeros((22, 2,3,2, len(bin_delay)))  
av_data_bin_resp  = np.zeros((22, 2,3,2, len(bin_resp)))  
av_data_bin_resp_lock = np.zeros((22, 2,3,2, len(bin_resp)))  


s = 0
for x in range(22):  
    nameSubj = 's' + str(x+1)
    
    # == Load individual data ==
    os.chdir(datapath_behav)
    data_behav = np.load(nameSubj +".npy") 

    coherence = data_behav[:, 5].astype(np.float64)
    stim_pres = data_behav[:, 6].astype(np.float64).astype(np.int64) 
    response_moment = data_behav[:,10].astype(np.float64).astype(np.int64)
    acc = data_behav[:, 11].astype(np.float64).astype(int) 
    rt = data_behav[:, 12].astype(np.float64)


    # == Load EMG artifact vector data by subject =
    vector_TT = np.load(datapath_vector_TT + "\s" + str(x+1) + "_vector_TT_EMG_EEG.npy")
        # (0 Response before Response Signal (RS))
        # (1 Missed Trial)
        # (2 Anticipation RT < 200 ms) 
        # (3 Could not detect EMG onsets)
        # (4 technical problem) 
        # (5 Valid trial, to include)
        # (6 Bad lag with corrected events)
        # (7 EEG eye artifact detected: blink, saccade)
        # (8 abnormal EEG voltage > 150µV)
    
    mask = (acc != -1) & (response_moment == 1) & (vector_TT == 5)
    
    for stim in range(2):
        for coh in range(3):
            for ch in range(2):
                temp_condition = data_bin_all[s, :, stim, coh, ch, :] # shape (960, x bin)
                mean_bin = temp_condition[mask & (stim_pres == possible_stim_levels[stim]) & (coherence == possible_coh_levels[coh])]
                # add mean of all trials in corresponding condition
                av_data_bin_all[s, stim, coh, ch, :] = np.nanmean(mean_bin, axis = 0)
                
                # -- idem with data per period ---
                temp_condition = data_bin_stim[s, :, stim, coh, ch, :] # shape (960, x bin)
                mean_bin = temp_condition[mask & (stim_pres == possible_stim_levels[stim]) & (coherence == possible_coh_levels[coh])]
                av_data_bin_stim[s, stim, coh, ch, :] = np.nanmean(mean_bin, axis = 0)
                
                temp_condition = data_bin_delay[s, :, stim, coh, ch, :] # shape (960, x bin)
                mean_bin = temp_condition[mask & (stim_pres == possible_stim_levels[stim]) & (coherence == possible_coh_levels[coh])]
                av_data_bin_delay[s, stim, coh, ch, :] = np.nanmean(mean_bin, axis = 0)
                
                temp_condition = data_bin_resp[s, :, stim, coh, ch, :] # shape (960, x bin)
                mean_bin = temp_condition[mask & (stim_pres == possible_stim_levels[stim]) & (coherence == possible_coh_levels[coh])]
                av_data_bin_resp[s, stim, coh, ch, :] = np.nanmean(mean_bin, axis = 0)
                
                temp_condition = data_bin_resp_lock[s, :, stim, coh, ch, :] # shape (960, x bin)
                mean_bin = temp_condition[mask & (stim_pres == possible_stim_levels[stim]) & (coherence == possible_coh_levels[coh])]
                av_data_bin_resp_lock[s, stim, coh, ch, :] = np.nanmean(mean_bin, axis = 0)

                assert mean_bin.shape[0] <= 160, f'More that 160 trials per conditions..., {mean_bin.shape[0]}'
    s+=1
print("Average completed")


# save result to plot in numpy
# os.chdir(datapath_JASP)
# np.save('av_data_bin_stim.npy', av_data_bin_stim)  #  22s, 2 stim, 3 coh, 2 channel emg, x bins
# np.save('av_data_bin_delay.npy', av_data_bin_delay) 
# np.save('av_data_bin_resp.npy', av_data_bin_resp) 
# np.save('av_data_bin_resp_lock.npy', av_data_bin_resp_lock) 




### ---------------------------------------------------- Convert and save long format to wide format for JASP analysis ----------------------------------------------------

def prepare_jasp_csv(df, bin_min=0, bin_max= None, channel=1, accuracy_exclude=-1, TT = 5):
    
    # ------------------ Filtrage ------------------
    df_filtered = df[
        (df['TT'] == TT) &
        (df['response_moment'] == 1) &
        (df['accuracy'] != accuracy_exclude) &
        (df['channel'] == channel) &
        (df['bin'] >= bin_min) &
        (df['bin'] <= bin_max)
    ].copy()
    
    # ------------------ Convertir en variables catégorielles ------------------
    df_filtered['subject'] = df_filtered['subject'].astype('category')
    df_filtered['trial'] = df_filtered['trial'].astype('category')

    # ------------------ Centrer les variables ------------------
    df_filtered['coherence_norm'] = df_filtered['coherence'] - df_filtered['coherence'].mean()
    df_filtered['stimulus_norm'] = df_filtered['stimulus'] - df_filtered['stimulus'].mean()
    df_filtered['PMA_norm'] = df_filtered['PMA']  # si pas de centrage requis

    # ------------------ Vérification des NAs et plage ------------------
    print(df_filtered.isna().sum())
    print(df_filtered['PMA'].min(), df_filtered['PMA'].max())

    # ------------------ Transformer PMA en variable ordinale ------------------
    df_filtered['PMA_ord'] = pd.Categorical(df_filtered['PMA'], categories=[0,1,2,3], ordered=True)

    # ------------------ Calculer les moyennes par condition pour RM-ANOVA ------------------
    df_means = df_filtered.groupby(
        ['subject', 'stimulus', 'coherence', 'bin'],
        as_index=False
    )['PMA_norm'].mean()

    # ------------------ Format wide pour JASP ------------------
    df_means['condition'] = (
        'stim' + df_means['stimulus'].astype(int).astype(str) +
        '_coh' + (df_means['coherence']*100).astype(int).astype(str) +
        '_bin' + df_means['bin'].astype(int).astype(str)
    )


    # ------------------ Moyenne par sujet × condition pour pivot ------------------
    df_means_subject = df_means.groupby(
        ['subject', 'condition'],
        as_index=False
    )['PMA_norm'].mean()
    

    df_wide = df_means_subject.pivot(index='subject', columns='condition', values='PMA_norm').reset_index()
    
    
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


# Select EMG channel ----
channel = 1 # 0 for opposite as response, 1 for same as resposne
if channel == 0: channel_name = '_oppCh'
else: channel_name = ''
os.chdir(datapath_JASP)

#  Stimulus period ------------------
df_wide = prepare_jasp_csv(df_long_stim, bin_min=0, bin_max=20, channel=channel)
df_wide.to_csv(f"PMA_rate_stim_JASP{channel_name}.csv", index=False)

#  Delay period ------------------
df_wide = prepare_jasp_csv(df_long_delay, bin_min=10, bin_max=30, channel=channel) # 1sec after stim offset (0-1s exclude see paper method)
df_wide.to_csv(f"PMA_rate_delay_JASP{channel_name}.csv", index=False)

# RS period ------------------
df_wide = prepare_jasp_csv(df_long_resp,bin_min=0, bin_max=4, channel=channel) # 1st after RS onset (+300ms)
df_wide.to_csv(f"PMA_rate_rs_JASP{channel_name}.csv", index=False)

# Response locked period ------------------
df_wide = prepare_jasp_csv(df_long_resp_onset,bin_min=0, bin_max=3, channel=channel) # 1st after RS onset (+300ms)
df_wide.to_csv(f"PMA_rate_resp_JASP{channel_name}.csv", index=False)



# [X]