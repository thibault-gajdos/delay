"""
EMG and EEG Experiment in Response signal variant random dot motion
- with 3 coherence levels and 2 durations of stimulus presentation (0.3s and 2s)
"""

# Importation -----
import os
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import arviz as az
import mne
import sys
from PyQt5.QtWidgets import QApplication, QDesktopWidget
import pandas as pd
import copy
current = os.path.join(os.getcwd()) 
toolbox_EMG_folder = r"C:\Users\erc\iCloudDrive\ACHARD_these\EMG_method\toolbox_EMG"
sys.path.extend([toolbox_EMG_folder]) # to add toolbox path to environment paths
import emgTools as emg
import evtTools as evt
from Events import Events
az.style.use('arviz-darkgrid')
np.set_printoptions(5, suppress=True)
os.chdir(current)


# Specify variable in function of your data --------------------------------
nb_subjects = 22

saving = False
possible_coh_levels = np.array([0, 0.12, 0.5 ])
possible_delay = np.array([18, 120]) # here delay is the stimulus lag presentation because delay fixed to 3s, in frame 18 for 300ms, 120 for 2s
nb_trial = 960
possible_level = 3*3 # 3mc x 3 stim pres

# RT trimming
trim_l = .200
nb_columns = 26 

### ------------------------------------------------------------------------------------------------------------

# Calculation of latencies from EMG (corrected_markers\segmented) or NPY (behavioral data)
latency = 'emg' # insert 'npy' or 'emg'

# Create empty containers to store data
mean_RTcorrect_per_coherence_delay = np.zeros((nb_subjects, possible_level))
mean_correct_per_coherence_delay = np.zeros((nb_subjects, possible_level))
proportion_correct_per_coherence_delay = np.zeros((nb_subjects, possible_level))
proportion_miss_per_coherence_delay = np.zeros((nb_subjects, possible_level))
proportion_beforeRS_per_coherence_delay = np.zeros((nb_subjects, possible_level)) #condition delay unqiuement ; réponses avant RS
proportion_anticipation_per_coherence_delay = np.zeros((nb_subjects, possible_level))

# New data for exportation
export = np.zeros((nb_trial*nb_subjects,nb_columns))

# data containers for latencies PRB
df_first_PRB_latencies = np.zeros((nb_trial*nb_subjects, 10))

# all latencies of the trial
df_latencies_PRB_same_chan_stim = np.zeros((nb_trial*nb_subjects, 4100))
df_latencies_PRB_opp_chan_PMT_stim = np.zeros((nb_trial*nb_subjects, 4100))
df_latencies_PRB_opp_chan_RT_stim = np.zeros((nb_trial*nb_subjects, 4100))

df_latencies_PRB_same_chan_delay = np.zeros((nb_trial*nb_subjects, 9100))
df_latencies_PRB_opp_chan_delay = np.zeros((nb_trial*nb_subjects, 9100))

df_latencies_PRB_same_chan_resp = np.zeros((nb_trial*nb_subjects, 10100))
df_latencies_PRB_opp_chan_PMT_resp = np.zeros((nb_trial*nb_subjects, 10100))
df_latencies_PRB_opp_chan_RT_resp = np.zeros((nb_trial*nb_subjects, 10100))


#------------ Data export parameters -------------
### Input yours folders
datapath_behav = os.path.abspath(os.path.join(current, "..", "data", "behavioral"))
datapath_EMG_signal = os.path.abspath(os.path.join(current, "..", "processing_emg", "data_saved")) # fif file step 1 emg preprocessing
datapath_EMG_segmented = os.path.abspath(os.path.join(current, "..", "data", "emg_segmented"))

# Channel 0 : left EMG / left
# Channel 1 : right EMG   / right


""" 

---- Columns output data .npy ----

|  | Columns Name           | Signification      | Var   | Details                                                                                                        |
|  |                        |                    |       |                                                                                                                |
|0 | subject_id             | ID Subject         | Int   | Subject number.                                                                                                |
|1 | trial_tot              | Trial number tot   | Int   | Number of trial in total for all block (1 -> 960 trials) 160 trials/conditions )                               |
|2 | block                  | Block number       | Int   | Number of trial in total for all block.                                                                        |
|3 | block_type             | Type of block      | Str   | There is 300 ms or 2s.                                                                                         |
|4 | trial_block            | Trial number block | Int   | Number of trial for block (1 -> 96 trials)                                                              |
|5 | condition_coherence    | Motion coherence   | Int   | Motion coherence of trial. Condition of coherence possible: 0, .12 or .50                                      |
|6 | stim_pres              | Time of pres. stim | Float | Time of presentation of stimulus during stimulus period: 300 ms or 2000ms                            |
|7 | direction_motion       | Dots direction     | Int   | Direction of dots in RDM for the trial. Left is 180, Right is 0.                                               |
|8 | key_required           | Expected key       | Int   | The key response expected for a correct response (1-left or 2-right).                                          |
|9 | key_pressed            | Pressedey          | Int   | The key response pressed (1-left or 2-right).                                                                  |
|10| response_moment        | Time of response   | Int   | The moment when the participant pressed a key can be too early (0), at the right time (1), or to late (-1).    |
|11| accuracy               | Accuracy           | Int   | The accuracy pf each subject for the task.                                                                     |
|12| rt                     | RT                 | Float | Response Time of participant from start of stimulus to button press (in ms), before the disapear of cue  .     |
|13| trigg_stim             | Stimulus trigger   | Int   | The stimulus trigger indicate conditions select of trial (delay, coherence, direction of dots). See * below.   |
|14| trigg_resp             | Response trigger   | Int   | The response trigger indicate response subject (correct or not, left or right). See ** below.                  |



---- New columns output data  ----

|--|-----------------------------------------------------------------------------|-------|-------------------------------------------------------------------------------------------------------------------------------------------|
|  | Columns Name                                                                | Var   | Details                                                                                                                                   |
|--|-----------------------------------------------------------------------------|-------|-------------------------------------------------------------------------------------------------------------------------------------------|
|0 | Subject_id                                                                  | Int   | Subject number, ID Subject.                                                                                                               |
|1 | Trial_id                                                                    | Int   | Number of trial in total for all block (1 -> 960 trials).                                                                                 |
|2 | Coherence                                                                   | Float | The three level of coherence: 0 / 0.12 / 0.50                                                                                             |
|3 | Accuracy                                                                    | Int   | 0 if incorrect, 1 if correct, -1 if no response.                                                                                          |
|4 | RT                                                                          | Float | Response time (sec). -1 if no response. RS correspond to response signal at repsonse.                                                     |
|5 | PMT                                                                         | Float | -1 if no response. PMT from RS to EMG onset of the response;                                                                              |
|6 | MT                                                                          | Float | -1 if no response. MT from EMG onset to response                                                                                          |
|7 | Trial_Type                                                                  | Int   | 1 missed trial / 2 anticipation / 3 could not detect EMG onsets / 4 technical problem / 5 valid trial / 6 bad lag bdf npy **See below.    |
|8 | Latency_first_partial_burst (stim onset -> resp)                            | Int   | -1 if no partial burst. The latency of first partial burst (between cue and response for DR).                                             |
|9 | Location_first_partial_burst (stim onset -> resp)                           | Int   | -1 if not EMG onsets/miss/anticipation trial, 0 if no partial burst, 1 if first partial on same side as the response, 2 if on opposite.   |
|10| motion direction                                                            | Int   | 180 if leftward, 0 if rightward.                                                                                                          |
|11| response_given                                                              | Int   | 1 if left response, 2 if right response.                                                                                                  |
|12| burst_between_EMGresp_and_resp                                              | Int   | 1 if partial burst between PMT and RT, 0 otherwise.                                                                                       |
|13| MAKES_PARTIAL_burst_between_EMGresp_and_resp                                | Int   | 1 if the partial burst between PMT and RT is the only partial burst in the trial.                                                         |
|14| stim pres (<=> delay presentation of stim)                                  | Int   | 300 or 2000 ms                                                                                                                            |
|15| response_moment                                                             | Int   | -1 if no response ; 0 if before Response Signal, 1 if after RS                                                                            |                                                                
|  |                                                                             |       |                                                                                                                                           | 
|16| nb_partial_burst_same_EMG_channel_as_response_during_stim_pres              | Int   | Number of partial burst in same channel as response during stimulus presentation (<2s). If FR condition, during PMT.                      |
|17| nb_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V1       | Int   | Number of partial burst in oppposite channel as response during stimulus presentation (<2s). If FR condition, during PMT.                 |
|18| nb_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V2       | Int   | Same as V1. You can change this to add FR condition, during RT to take into account synesis) [here RS blocks, V1 and V2 are identical]    |
|19| nb_partial_burst_same_EMG_channel_during_delay                              | Int   | Number of partial burst in correct channel during delay period                                                                            |
|20| nb_partial_burst_opposite_EMG_channel_during_delay                          | Int   | Number of partial burst in incorrect channel during delay period.                                                                         |
|21| nb_partial_burst_same_EMG_channel_as_response_during_response_signal        | Int   | Number of partial burst in correct channel during response signal period during PMT.                                                      |
|22| nb_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V1 | Int   | Number of partial burst in incorrect channel during response signal period during PMT.                                                    |
|23| nb_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V2 | Int   | Number of partial burst in incorrect channel during response signal period during RT to take into account syncinesis.                     |
|24| Latency_first_partial_burst_RS  (RS -> resp)
|25| Location_first_partial_burst_RS (RS -> resp)
|--|-----------------------------------------------------------------------------|-------|-------------------------------------------------------------------------------------------------------------------------------------------|



* Stimulus Onset Trigger for 3s delay
Trigger Stim code:

|                 | Stim Pres |  coherence  0.  |  coherence .12  |  coherence .50  |  
|                 |           |                 |                 |                 |  
| RS 300          |  300 ms   |   L 21 / R 24   |   L 22 / R 25   |   L 23 / R 26   |  
| RS 2000         |  2000 ms  |   L 31 / R 34   |   L 32 / R 35   |   L 33 / R 36   |  

** Response Signal trigger (after disapear of cross = trigg_stim + 100)



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





** Trial_Type:
(0 Response before Response Signal (RS))
(1 Missed Trial)
(2 Anticipation RT < 200 ms) 
(3 Could not detect EMG onsets)
(4 technical problem, remove this trial) #will override all other trial type
(5 Valid trial, no problem cited above)
(6 Bad lag with corrected events)


Note: 
left chan = 0
right chan = 1

"""

coh = possible_coh_levels.shape[0]
stim_lag = possible_delay.shape[0]

ind_av_EMG_Correct =np.zeros((stim_lag, coh, nb_subjects, 412))
ind_av_EMG_Correct_resp =np.zeros((stim_lag, coh,  nb_subjects, 361))

ind_av_EMG_All =np.zeros((stim_lag, coh, nb_subjects, 412))
ind_av_EMG_All_resp =np.zeros((stim_lag, coh,  nb_subjects, 361))

data = {
    'subject': [0],
    'trial': [0],
    'coherence': [0],
    'lag':[0],
    'period': [0],
    'chan_same_opp_resp': [0],
    'latence': [0],
    'PMA': [0],
    '1st PRB': [0],
    'trial_type': [0], 
    'accuracy': [0]
}
df_PRB_analysis_exp4 = pd.DataFrame(data)



# Navigate between subjects --------------------------------------------------------------------
count = 0
for s in range(nb_subjects):
    nameSubj = 's' + str(s+1)

    # Load individual data ----------------------------------------------
    data_behav = np.load(os.path.join(datapath_behav, nameSubj +".npy")) 

    coherence = data_behav[:, 5].astype(np.float64)
    stim_pres = data_behav[:, 6].astype(np.float64).astype(np.int64) 
    response_moment = data_behav[:,10].astype(np.float64).astype(np.int64)
    acc = data_behav[:, 11].astype(np.float64).astype(int) 
    rt = data_behav[:, 12].astype(np.float64)
    motion_direction = data_behav[:, 7].astype(np.float64).astype(np.int64)  # 180: left   0: right 
    key_pressed = data_behav[:,9].astype(np.float64).astype(np.int64) 
    trigg_resp = data_behav[:, 14].astype(np.float64).astype(int)
    
    
    # Here delay variable is the delay of presentation of the stimulus during stimulus presentation period (also called stimLag)
    delay = np.zeros((len(data_behav)))
    for i in range(len(data_behav)):
        if data_behav[i][6] == '18': delay[i] = 0.300
        elif data_behav[i][6] == '120': delay[i] = 2.000
    delay = delay.astype(np.float64)

    
    #Load EMG data -----------------------------------------------------------
    data_fname = os.path.join(datapath_EMG_signal, nameSubj + '_EMGstep1.fif')
    raw = mne.io.Raw(data_fname, preload=True)
    sf = raw.info['sfreq']
    EMG_data = raw.get_data(picks=['EMGL','EMGR']) # filtered
    EMG_data = EMG_data * 1000000 
    EMG_data = np.absolute(EMG_data) # rectify EMG_data
    
    #Load segmented events ----------------------------------------------------
    sf = 1024.0
    MARK_seg = evt.loadSegmented(os.path.join(datapath_EMG_segmented, nameSubj + "_segmentedEvts.csv"),headerlines=1, sep=',', sf=sf,
                        colSample=0, colTime=1, colCode=2, colChan=3,
                        colTrialIdx=4, 
                        colTminSample=5, colTminTime=6,
                        colT0Sample=7, colT0Time=8, 
                        colTmaxSample=9, colTmaxTime=10)
    MARK_seg_modified = copy.deepcopy(MARK_seg)


#     NOTES:
#         MARK_seg.listEvtsTrials  epEvts.listEvtsTrials[e] contains triggers and corresponding latencies for epoch e
#         MARK_seg.listEvtsTrials[0].code  returns np.array containg triggers of epoch e
#         MARK_seg.listEvtsTrials[0].Latency.sample  returns np.array containg latency (in samples) of triggers of epoch e
#         MARK_seg.listEvtsTrials[0].Latency.time  returns np.array containg latency (in ms) of triggers of epoch e
#         MARK_seg.listEvtsTrials[0].findEvents(sample = OR  time = OR code= OR channel)Finds the specified event(s) and returns the position(s).
#         MARK_seg.nbTrials()   interger number of epochs

#         >> Relationship with continuous signal: 
#         MARK_seg.tmin.time np vector array   tmin[n] contains starting latency of epoch n in continuous signal
#         MARK_seg.t0.time np vector array   t0[n] contains latency of locking event of epoch n in continuous signal
#         MARK_seg.tmax.time np vector array   tmax[n] contains ending latency of epoch n in continuous signal


    """ 
    Trigger Stim code:
    |                 | Stim Pres |  coherence  0.  |  coherence .12  |  coherence .50  |  
    |                 |           |                 |                 |                 |  
    | RS 300          |  300 ms   |   L 21 / R 24   |   L 22 / R 25   |   L 23 / R 26   |  
    | RS 2000         |  2000 ms  |   L 31 / R 34   |   L 32 / R 35   |   L 33 / R 36   |  

    """
    
    # ID triggers
    emg_channels_idx = {0: 'EMGL', 1: 'EMGR'}
    trigg_stim_offset = 1
    trigg_stim_offset_id = {'stim_offset':1}
    
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

    trigg_resp_id = {'left_correct': 101, 'left_incorrect': 100, 'right_correct' : 103, 'right_incorrect':102, 
                     'left_correct_beforeRS': 201, 'left_incorrect_beforeRS': 200, 'right_correct_beforeRS' : 203, 'right_incorrect_beforeRS':202, 
                     'no_resp':250}
    
    trigg_emg_id = {'onset': 350 , 'offset': 351 }



    # ----------------------------------------- For each trial (i) per subject (s) -----------------------------------------------------------
    for i in range(MARK_seg.nbTrials()):

        # Count number of EMG onsets in the current trial (before the response), regardless of the channel
        len_EMG_onsets = len(MARK_seg.listEvtsTrials[i].code[MARK_seg.listEvtsTrials[i].code == 350])

        # For all response motion and key pressed
        export[count,10] = motion_direction[i]
        export[count,11] = key_pressed[i]
        export[count,14] = delay[i]
        export[count,7] = -1

        # for all trial, can be change below
        export[count,24] = -1 # no 1st letency PRB
        export[count,25] = -1 # no location to 1st PRB 


        # Depending on the type of trigger response: ------------------------------------------
        # Missed trial 
        if trigg_resp[i] == 250: 
            export[count,0]= s+1
            export[count,1]= i + 1 # trial
            export[count,2] = coherence[i]
            export[count,3:7] = -1 # Acc, RT, PMT, MT
            export[count,7] = 1 # missed trial
            export[count,8] = -1 # latency first burst
            export[count,9] = -1 # location first burst
            export[count,12] = -1 # burst_between_EMGresp_and_resp 
            export[count,13] = -1 # MAKES_PARTIAL_burst_between_EMGresp_and_resp 
            export[count,15] = -1 #response moment
            assert response_moment[i] == -1, print ("response_moment not congruent with trigger 250")
  
        # All trials with response --------------------------------------------------------
        if (trigg_resp[i] != 250):
            #Could NOT detect EMG onsets
            if (len_EMG_onsets == 0):
                export[count,0]= s+1
                export[count,1]= i + 1
                export[count,2] = coherence[i]
                export[count,3] = acc[i]
                export[count,4] = rt[i]
                export[count,5] = -1
                export[count,6] = -1
                export[count,7] =  3 #could not detect EMG onsets            
                export[count,8] = -1 #latency first burst
                export[count,9] = -1 #location first burst  
                export[count,12] = -1# burst_between_EMGresp_and_resp 
                export[count,13] = -1# MAKES_PARTIAL_burst_between_EMGresp_and_resp 
                export[count,15] = response_moment[i]              
                     
            #Anticipations -> if want anticipation
            if (rt[i] < trim_l):
                export[count,0]= s+1
                export[count,1]= i + 1
                export[count,2] = coherence[i]
                export[count,3:7] = -1 #Acc, RT, PMT, MT
                export[count,7] = 2 #anticipation         ####### will override category 3   
                export[count,8] = -1 #latency first burst
                export[count,9] = -1 #location first burst  
                export[count,12] = -1# burst_between_EMGresp_and_resp 
                export[count,13] = -1# MAKES_PARTIAL_burst_between_EMGresp_and_resp 
                export[count,15] = response_moment[i]  
            

            # Detection EMG onsets and response_moment ok ------------------------------------------
            if (len_EMG_onsets > 0) and (response_moment[i]==1) and (rt[i] >= trim_l): 

                #### Number of partial during stimulus presentation period (egal to len all burst below)
                # latency partial burst during stimulus presentation period
                lat_partial_burst_same_EMG_channel_as_response_during_stim_pres = []      
                lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V1 = [] 
                lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V2 = [] #[dans DR blocks, V1 et V2 sont identiques]

                # latency partial burst during delay period
                lat_partial_burst_same_EMG_channel_during_delay  = []  
                lat_partial_burst_opposite_EMG_channel_during_delay = [] 

                # latency partial burst during response signal period
                lat_partial_burst_same_EMG_channel_as_response_during_response_signal = []      #(pendant PMT
                lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V1 = [] #(pendant PMT
                lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V2 = [] #(pendant RT pour prendre en compte syncinésies


                # ---------------------------------------------------------------------------------------------------------------#
                ########### -------------------------- Trigger response (100,101,102,103) ---------------------------- ###########
                # ---------------------------------------------------------------------------------------------------------------#
                
                ### Complete export data in function of type response for trial with trigger response  
                """ Memo:
                Incorrect left  response 100  / Correct left  response 101
                Incorrect right response 102  / Correct right response 103
                """ 
    
                if (trigg_resp[i] in [100, 101, 102, 103, 99]): # To be sure that we are in the coded answers    
                
                    export[count,0]= s+1
                    export[count,1]= i + 1
                    export[count,2] = coherence[i]

                    if trigg_resp[i] in [100, 102]:   export[count,3] = 0 #Acc
                    elif trigg_resp[i] in [101, 103]: export[count,3] = 1
                    else:  export[count,3] = -1

                    export[count,7] = 5 # Pour les essais valides, car sinon la valeur de base est 0
                    export[count,15] = response_moment[i]   
    
                    # idx stim is RS signal stim here 
                    idx_stim = int(MARK_seg.listEvtsTrials[i].findEvents(code=list(trigg_RS_id.values()))) # second trigg_stim_signal_response after RS (clock was reset at RS)
                    idx_initial_stim = int(MARK_seg.listEvtsTrials[i].findEvents(code=list(trigg_stim_id.values()))) # first stim onset

                    # Find latency of stim RS (and response after)
                    latency_stim = MARK_seg.listEvtsTrials[i].Latency.time[idx_stim]


                    # if technical problems
                    if (s==99): 
                        print('add technical problems')
                    elif (MARK_seg.listEvtsTrials[i].code[-2] == 99) or (MARK_seg.listEvtsTrials[i].code[-4] == 99): # trial with lag to exclude
                        export[count,0]= s+1
                        export[count,1]= i + 1 # trial
                        export[count,2] = coherence[i]
                        export[count,3:7] = -1 # Acc, RT, PMT, MT
                        export[count,7] = 6   # bad lag trial (compute on correct events)
                        export[count,8] = -1 # latency first burst
                        export[count,9] = -1 # location first burst
                        export[count,12] = -1 # burst_between_EMGresp_and_resp 
                        export[count,13] = -1 # MAKES_PARTIAL_burst_between_EMGresp_and_resp 
                        export[count,15] = -1 #response moment
                        print('Bad lag continue...')
                        print(MARK_seg.listEvtsTrials[i])
                        count = count + 1 #trial count  
                        continue
                    else:                 
                        idx_resp = int(MARK_seg.listEvtsTrials[i].findEvents(code=trigg_resp[i])) # Index of trigger response
                        # check RS trigger
                        if (response_moment[i] == 1):
                            for rs in range(len(trigg_RS)):
                                if trigg_RS[rs] in MARK_seg.listEvtsTrials[i].code:
                                    trial_trigg_RS = trigg_RS[rs]
                            idx_RS = int(MARK_seg.listEvtsTrials[i].findEvents(code=trial_trigg_RS)) # Index
                            assert idx_RS == idx_stim, print('RS index not congruent')
                            
                        else:
                            idx_RS = 0 # because not RS in this case
                            
                    nb_trigg= np.flip(np.arange(idx_resp)) # reverse trigger to check EMG onset before response. So first trigger is always 350 (the two next are response and 351)

                    # Find Latency
                    if latency == 'emg':
                        lat_resp = MARK_seg.listEvtsTrials[i].Latency.time[idx_resp] - MARK_seg.listEvtsTrials[i].Latency.time[idx_stim] # EMG latency of response : time response - time RS
                        export[count,4] = lat_resp #RT
                    elif latency == 'npy':
                        lat_resp = rt[i] #MARK_seg.listEvtsTrials[i].Latency.time[idx_resp] # behavior latency
                    
    
                    # --- Find last EMG onset in the incorrect channel ---
                    # channel_associated_with_response is the channel of the response
                    if trigg_resp[i] in [100,101]:
                        channel_associated_with_response = 0 # for incorrect/correct left response channel of the response is left (0)
                        channel_not_associated_with_response = 1
                    elif trigg_resp[i] in [102,103]:
                        channel_associated_with_response = 1  # for incorrect/correct left response channel of the response is right (1)
                        channel_not_associated_with_response = 0    
                    check = 0
                    counter = 0
                    
                    # --- Find first EMG onset correspond to EMG onset of the response and break ---
                    for n in nb_trigg:
                        if  MARK_seg.listEvtsTrials[i].code[n] == 350:
                            if (MARK_seg.listEvtsTrials[i].chan[n] == channel_associated_with_response): 
                                check = 1
                                idx_pmt_burst = n # save position of last burst of PMT
                                lat_burst_LAST = MARK_seg.listEvtsTrials[i].Latency.time[n]

                                export[count,5] = lat_burst_LAST - latency_stim #PMT   -note that latency_stim != 0 here
                                export[count,6] = lat_resp - (lat_burst_LAST - latency_stim) # MT  

                                if counter > 0:
                                    export[count,12] = 1 #During the EMG of the response, there is an activity in the opposite response channel
                                else:
                                    export[count,12] = 0  

                                                                            
                                # RECODE EMG onset marker ----------------------------------
                                if trigg_resp[i] == 100:#Incorrect left  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1300
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1301
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1302  
                                          
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1700
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1701
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1702   

                                if trigg_resp[i] == 101:#Correct left  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1300
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1301
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1302  
                                             
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1700
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1701
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1702   

                                if trigg_resp[i] == 102:#Incorrect right  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3300
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3301
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3302  
                                             
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3700
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3701
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3702  

                                if trigg_resp[i] == 103:#Correct right  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3300
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3301
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3302  
                                             
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3700
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3701
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3702                                                                                          
                                break
                            
                        counter = counter + 1                               
                    if (check == 0) and (export[count,7] != 4): # if trial  is already classified as a technical problem, it will not be taken into account
                        print ("WARNING: failed detecting last EMG onset leading to response - Trial:",i)     
                        print("MARK_seg.listEvtsTrials[i]:",MARK_seg.listEvtsTrials[i]) 
                        
                        
                    # Now transform response trigger ----------------------------------------------
                    if trigg_resp[i] == 100:#Incorrect left  response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100002
                                

                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120002                                       
                                
                    if trigg_resp[i] == 101:#Correct left  response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100002
                                

                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120002
                                
                                
                    if trigg_resp[i] == 102:#Incorrect right response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160002
                                
                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180002                                    
                                
                    if trigg_resp[i] == 103:#Correct right response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160002
                                
                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180000
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180001
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180002   

                        # MARK_seg.listEvtsTrials[i].code          = array([ 31,   1, 131, 350, 101, 351])
                        # MARK_seg_modified.listEvtsTrials[i].code = array([    31,      1,    131,   2300, 130000,    351])


                    
                    
                    
                    # ---- Mean 1st PRB RS to response ----
                    if (response_moment[i] == 1):
                        # -- in response period only 
                        first_partial_burst_location = 0
                        for b in range(idx_RS, idx_pmt_burst):#idx_resp corresponds to index of response and -1 to not count the 350 onset of response
                            time = MARK_seg.listEvtsTrials[i].Latency.time[b]

                            # Location and lantecy of first burst all period
                            if (MARK_seg.listEvtsTrials[i].code[b] == 350):
                                if first_partial_burst_location == 0:
                                    export[count,8] = MARK_seg.listEvtsTrials[i].Latency.time[b] - MARK_seg.listEvtsTrials[i].Latency.time[idx_stim]  # correspond to time withou stim and delay period
                                    first_partial_burst_location = 1
                                    
                                    # Same channel
                                    if MARK_seg.listEvtsTrials[i].chan[b] == channel_associated_with_response:
                                        export[count,9] = 1 #first partial burst detected on the same side as the response  
                                    # Opposite channel
                                    if MARK_seg.listEvtsTrials[i].chan[b] == channel_not_associated_with_response:
                                        export[count,9] = 2 #first partial burst detected on the opposite side as the response
                                        

                        #--  Mean 1st PRB all trial stim-resp ----
                        first_partial_burst_location = 0
                        for b in range(idx_resp-2):#idx_resp corresponds to index of response and -1 to not count the 350 onset of response
                            time = MARK_seg.listEvtsTrials[i].Latency.time[b]

                            # Location and lantecy of first burst
                            if (MARK_seg.listEvtsTrials[i].code[b] == 350):
                                if first_partial_burst_location == 0:
                                    export[count,24] = MARK_seg.listEvtsTrials[i].Latency.time[b] 
                                    first_partial_burst_location = 1
                                    # Same channel
                                    if MARK_seg.listEvtsTrials[i].chan[b] == channel_associated_with_response:
                                        export[count,25] = 1 #first partial burst detected on the same side as the response  
                                    # Opposite channel
                                    if MARK_seg.listEvtsTrials[i].chan[b] == channel_not_associated_with_response:
                                        export[count,25] = 2 #first partial burst detected on the opposite side as the response


                    # Find stimulus offset index
                    idx_offset = int(MARK_seg.listEvtsTrials[i].findEvents(code=trigg_stim_offset))  # must be equal to int(1)
                                      

                    for b in range(idx_stim): #idx_stim corresponds to index of RS trigger
                        time = MARK_seg.listEvtsTrials[i].Latency.time[b] # e.g. latency of RS 

                        # For stimulus presentation period
                        # if time < export[count,14]: # if time < stim duration (.3 /  2.)
                        if (b <= idx_offset):
                            if (MARK_seg.listEvtsTrials[i].code[b] == 350):
                                if MARK_seg.listEvtsTrials[i].chan[b] == channel_associated_with_response:
                                    lat_partial_burst_same_EMG_channel_as_response_during_stim_pres.append(MARK_seg.listEvtsTrials[i].Latency.time[b])
                                if MARK_seg.listEvtsTrials[i].chan[b] == channel_not_associated_with_response:
                                    lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V1.append(MARK_seg.listEvtsTrials[i].Latency.time[b])
                                    lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V2.append(MARK_seg.listEvtsTrials[i].Latency.time[b])
                        
                        # For delay period
                        # if (time >= export[count,14]) and (time < (export[count,14] + 3)): #export[count,14] is stim duration + 3s delay
                        elif (b > idx_offset):    
                            if (MARK_seg.listEvtsTrials[i].code[b] == 350):
                                if MARK_seg.listEvtsTrials[i].chan[b] == channel_associated_with_response:
                                    lat_partial_burst_same_EMG_channel_during_delay.append(MARK_seg.listEvtsTrials[i].Latency.time[b])
                                    # debug
                                    # if (MARK_seg.listEvtsTrials[i].Latency.time[b] < 0):
                                    #     print(MARK_seg.listEvtsTrials[i])
                                        
                                if MARK_seg.listEvtsTrials[i].chan[b] == channel_not_associated_with_response:
                                    lat_partial_burst_opposite_EMG_channel_during_delay.append(MARK_seg.listEvtsTrials[i].Latency.time[b])

                    # For response period during PMT (not take account 350 onset of response)
                    for r in range(idx_stim, idx_pmt_burst):  # from RS trigger to EMG onset of response (response not take account)
                        if (MARK_seg.listEvtsTrials[i].code[r] == 350):
                                if MARK_seg.listEvtsTrials[i].chan[r] == channel_associated_with_response:
                                    lat_partial_burst_same_EMG_channel_as_response_during_response_signal.append(MARK_seg.listEvtsTrials[i].Latency.time[r])
                                if MARK_seg.listEvtsTrials[i].chan[r] == channel_not_associated_with_response:
                                    lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V1.append(MARK_seg.listEvtsTrials[i].Latency.time[r])
                                        
                    # For response period during RT # check in oppsite channel as response if there is a PMA
                    for r in range(idx_stim, idx_resp):
                        if (MARK_seg.listEvtsTrials[i].code[r] == 350):
                            if MARK_seg.listEvtsTrials[i].chan[r] == channel_not_associated_with_response:
                                lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V2.append(MARK_seg.listEvtsTrials[i].Latency.time[r])

                    # ---- Save results as number per perriod ----
                    # Burst during stim presentation
                    export[count,16]= len(lat_partial_burst_same_EMG_channel_as_response_during_stim_pres)
                    export[count,17]= len(lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V1) # V1 and V2 is identical here, because only RS condition 
                    export[count,18]= len(lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V2)

                    # Burst during delay period
                    export[count,19]= len(lat_partial_burst_same_EMG_channel_during_delay)
                    export[count,20]= len(lat_partial_burst_opposite_EMG_channel_during_delay)

                    #  Burst during response signal
                    export[count,21]= len(lat_partial_burst_same_EMG_channel_as_response_during_response_signal)      # during PMT

                    export[count,22]= len(lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V1) # during PMT
                    export[count,23]= len(lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V2) # during RT (to take into account PMA in oppposite channel as response)

                    # export 13. For RS, we check if there's only one onset 350 in the response period (instead of making a list, I've added in MARK ....code[idx_stim:idx_resp] 
                    # to be in the right period. It works well, you can test with s = 8 which ends with a block RS
                    if list(MARK_seg.listEvtsTrials[i].code[idx_stim:idx_resp] == 350).count(True) == 1:
                        export[count,13] = 1 #the partial burst between PMT and RT is the only partial burst in the trial.          

            
            
            
            
            
            
            
            
            # Detection EMG onsets and response_moment before RS ------------------------------------------------------------------
            if (len_EMG_onsets > 0) and (response_moment[i]==0):
                # ---------------------------------------------------------------------------------------------------------------#
                ########### ---------------------------- Premature response (200,201,202,203) ------------------------ ###########
                # ---------------------------------------------------------------------------------------------------------------#
                
                if (trigg_resp[i] in [200,201,202,203]): # To be sure that we are in the coded answers    
                
                    export[count,0]= s+1
                    export[count,1]= i + 1
                    export[count,2] = coherence[i]

                    if trigg_resp[i] in [200, 202]:   export[count,3] = 0 #Acc
                    elif trigg_resp[i] in [201, 203]: export[count,3] = 1
                    else:  export[count,3] = -1

                    export[count,7] = 0 # Pour les essais valides premature response
                    export[count,15] = response_moment[i]   
    
                    # idx stim is first stimulus for premature response
                    idx_stim = int(MARK_seg.listEvtsTrials[i].findEvents(code=list(trigg_stim_id.values()))) # first stim onset

                    # Find latency of stim and response & Index of the response
                    latency_stim = MARK_seg.listEvtsTrials[i].Latency.time[idx_stim]


                    # if technical problems
                    if (s==99): print('add technical problems')
                    else:                 
                        idx_resp = int(MARK_seg.listEvtsTrials[i].findEvents(code=trigg_resp[i])) # Index of trigger response
                        idx_RS = 0
                    nb_trigg= np.flip(np.arange(idx_resp))
    

                    if latency == 'emg':
                        lat_resp = MARK_seg.listEvtsTrials[i].Latency.time[idx_resp] - MARK_seg.listEvtsTrials[i].Latency.time[idx_stim] # EMG latency 
                        export[count,4] = lat_resp #RT
                    elif latency == 'npy':
                        lat_resp = rt[i] #MARK_seg.listEvtsTrials[i].Latency.time[idx_resp] # behavior latency
                    
    
                    # --- Find last EMG onset in the incorrect channel ---

                    # channel_associated_with_response is the channel of the response
                    if trigg_resp[i] in [200,201]:
                        channel_associated_with_response = 0 # for incorrect/correct left response channel of the response is left (0)
                        channel_not_associated_with_response = 1
                    elif trigg_resp[i] in [202,203]:
                        channel_associated_with_response = 1  # for incorrect/correct left response channel of the response is right (1)
                        channel_not_associated_with_response = 0    
                    check = 0
                    counter = 0
                    

                    for n in nb_trigg:
                        if  MARK_seg.listEvtsTrials[i].code[n] == 350:
                            if (MARK_seg.listEvtsTrials[i].chan[n] == channel_associated_with_response): 
                                check = 1
                                idx_pmt_burst = n # save position of last burst of PMT
                                lat_burst_LAST = MARK_seg.listEvtsTrials[i].Latency.time[n]

                                export[count,5] = lat_burst_LAST - latency_stim #PMT   -note that latency_stim != 0 here
                                export[count,6] = lat_resp - (lat_burst_LAST - latency_stim) # MT  

                                if counter > 0:
                                    export[count,12] = 1 #During the EMG of the response, there is an activity in the opposite response channel
                                else:
                                    export[count,12] = 0  
                                    
                                          
                                #RECODE EMG onset marker -----------------------------------------------
                                if trigg_resp[i] == 200:#Incorrect left  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1300 + 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1301+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1302  + 100
                                       
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1700+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1701+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 1702   + 100

                                if trigg_resp[i] == 201:#Correct left  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2300+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2301+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2302  + 100
                                         
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2700+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2701+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 2702   + 100

                                if trigg_resp[i] == 202:#Incorrect right  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3300+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3301+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3302  + 100
                                         
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3700+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3701+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 3702  + 100

                                if trigg_resp[i] == 203:#Correct right  response
                                    if export[count,14] == .300:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4300+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4301+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4302  + 100
                                        
                                    if export[count,14] == 2.:
                                        if coherence[i] == 0:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4700+ 100
                                        if coherence[i] == 0.12:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4701+ 100
                                        if coherence[i] == 0.5:
                                            MARK_seg_modified.listEvtsTrials[i].code[n] = 4702   + 100                                                                                       
                                break
                                    
                        counter = counter + 1 
                                                      
                    
                    if (check == 0) and (export[count,7] != 4): # if trial  is already classified as a technical problem, it will not be taken into account
                        print ("WARNING: failed detecting last EMG onset leading to response - Trial:",i)     
                        print("MARK_seg.listEvtsTrials[i]:",MARK_seg.listEvtsTrials[i]) 
                    
                    
                    #now transform response trigger ----------------------------------------------------
                    if trigg_resp[i] == 200:#Incorrect left  response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 100002+ 100
                                
                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 120002   + 100                                    
                                
                    if trigg_resp[i] == 201:#Correct left  response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 130000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 130001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 130002+ 100
                        
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 150000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 150001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 150002+ 100
                                
                                
                    if trigg_resp[i] == 202:#Incorrect right response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 160002+ 100
                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 180002   + 100                                 
                                
                    if trigg_resp[i] == 203:#Correct right response
                        if export[count,14] == .300:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 190000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 190001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 190002+ 100
                                
                        if export[count,14] == 2.:
                            if coherence[i] == 0:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 210000+ 100
                            if coherence[i] == 0.12:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 210001+ 100
                            if coherence[i] == 0.5:
                                MARK_seg_modified.listEvtsTrials[i].code[idx_resp] = 210002   + 100
        
        
        # --- simply name of variable ---
        #stim
        stim_same = lat_partial_burst_same_EMG_channel_as_response_during_stim_pres        # is latence of detected partial EMG
        stim_opp = lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V1  # V1 and V2 is identical
        stim_opp_bis = lat_partial_burst_opposite_EMG_channel_as_response_during_stim_pres_V2

        #delay
        delay_same = lat_partial_burst_same_EMG_channel_during_delay
        delay_opp = lat_partial_burst_opposite_EMG_channel_during_delay

        #resp
        resp_same = lat_partial_burst_same_EMG_channel_as_response_during_response_signal
        resp_opp = lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V1
        resp_opp_bis = lat_partial_burst_opposite_EMG_channel_as_response_during_response_signal_V2


        # ------ Data container for PRB analysis ------
        period_name = ['stim', 'delay', 'resp']
        ch_names = ['same', 'opp']
        select_data_lat = [[stim_same, stim_opp], [delay_same, delay_opp], [resp_same, resp_opp]]

        for period in range(len(period_name)):
            for ch in range(len(ch_names)):

                data = select_data_lat[period][ch]

                # if PRB, add info in df_PRB_analysis 'subject  trial  period  chan_same_opp_resp  latence'
                if len(data) == 1 :
                    period_n = [period_name[period]]
                    chan = [ch_names[ch]]
                    trial_type = export[count,7]
                    acc_trial = export[count,3]
                    
                    if period_n == 'resp': first_PRB = 1
                    else: first_PRB = 0

                    dataPRB_trial = {'subject': s+1,'trial': i+1,'coherence': coherence[i],'lag': delay[i],'period': period_n,'chan_same_opp_resp': chan,
                                     'latence': data[0], 'PMA': 1, '1st PRB': first_PRB, 'trial_type': trial_type, 'accuracy': acc_trial}
                    new_row = pd.DataFrame(dataPRB_trial)
                    df_PRB_analysis_exp4 = pd.concat([df_PRB_analysis_exp4, new_row], ignore_index=True)

                elif len(data) > 1 :                  
                    subj_list = np.full(len(data), s+1) # creer une liste de plusieurs fois le numero de sujet
                    trial_list = np.full(len(data), i+1) 
                    coh_list = np.full(len(data), coherence[i]) 
                    delay_list = np.full(len(data), delay[i])   # is the lag
                    period_list = np.full(len(data), period_name[period])
                    chan_list = np.full(len(data), ch_names[ch])
                    lat_list = np.full(len(data), data)
                    n_PMA_list = np.full(len(data), 1)
                    
                    if period_list[0] == 'resp': prb_1st_list = [1]
                    else: prb_1st_list = [0]
                    prb_1st_list.extend(np.full(len(data)-1, 0))
                    
                    trial_type_list = np.full(len(data), export[count,7])
                    acc_list = np.full(len(data), export[count,3])

                    dataPRB_trial = {'subject': subj_list,'trial': trial_list,'coherence': coh_list,'lag': delay_list, 'period': period_list,
                                     'chan_same_opp_resp': chan_list,'latence': lat_list, 'PMA': n_PMA_list, '1st PRB': prb_1st_list, 'trial_type': trial_type_list, 'accuracy': acc_list}
                    new_rows = pd.DataFrame(dataPRB_trial)
                    df_PRB_analysis_exp4 = pd.concat([df_PRB_analysis_exp4, new_rows], ignore_index=True)
                    
                else:
                    period_n = [period_name[period]]
                    chan = [ch_names[ch]]
                    trial_type = export[count,7]
                    acc_trial = export[count,3]
                    
                    dataPRB_trial = {'subject': s+1,'trial': i+1,'coherence': coherence[i],'lag': delay[i],'period': period_n,'chan_same_opp_resp': chan,
                                     'latence': -1, 'PMA': 0, '1st PRB': -1, 'trial_type': trial_type, 'accuracy': acc_trial}
                    new_row = pd.DataFrame(dataPRB_trial)
                    df_PRB_analysis_exp4 = pd.concat([df_PRB_analysis_exp4, new_row], ignore_index=True)


        count = count + 1 #trial count                                      

    MARK_seg_modified_continuous = MARK_seg_modified.asContinuous()[0] 

    # Now segment trials locked on EMG data -----------------------------------------------------------
    locking_events_lag300 = np.array([ [1300,3300],
                                               [1301,3301],
                                               [1302,3302]])

    locking_events_lag2000 = np.array([ [1700,3700],
                                               [1701,3701],
                                               [1702,3702]])


    resp_locking_events_lag300 = np.array([ [100000,160000],
                                               [100001,160001],
                                               [100002,160002]])  

    resp_locking_events_lag2000 = np.array([ [120000,180000],
                                               [120001,180001],
                                               [120002,180002]])    

    
    # Parameters --------------------
    min_stim = -.1
    max_stim = .3
    epochTime = evt.times(min_stim, max_stim, sf)   
    t0 = emg.findTimes(0, epochTime)    
    
    min_resp = -.3
    max_resp = .05
    epochTime_resp = evt.times(min_resp, max_resp, sf)      
    


    for i in range (3): # coherence
        
        # --- For incorrect + correct trials ---
        # 300ms stimulus duration 
        # EMG-locked -----------
        Epochs_L = MARK_seg_modified_continuous.segment(codeT0=locking_events_lag300[i,0], tmin=min_stim, tmax=max_stim)
        EMG_Epochs_L = evt.getDataEpochs(EMG_data, Epochs_L)#(48, 2, 412)
        EMG_Epochs_L =EMG_Epochs_L[:,0,:]    
         
        Epochs_R = MARK_seg_modified_continuous.segment(codeT0=locking_events_lag300[i,1], tmin=min_stim, tmax=max_stim)
        EMG_Epochs_R = evt.getDataEpochs(EMG_data, Epochs_R)    
        EMG_Epochs_R =EMG_Epochs_R[:,1,:]    

        EMG_Epochs_av = np.concatenate((EMG_Epochs_L,EMG_Epochs_R)).mean(axis = 0)        
    
        # response-locked -----------
        Epochs_L_resp = MARK_seg_modified_continuous.segment(codeT0=resp_locking_events_lag300[i,0], tmin=min_resp, tmax=max_resp)
        EMG_Epochs_L_resp = evt.getDataEpochs(EMG_data, Epochs_L_resp)
        EMG_Epochs_L_resp =EMG_Epochs_L_resp[:,0,:]
    
        Epochs_R_resp = MARK_seg_modified_continuous.segment(codeT0=resp_locking_events_lag300[i,1], tmin=min_resp, tmax=max_resp)
        EMG_Epochs_R_resp = evt.getDataEpochs(EMG_data, Epochs_R_resp)    
        EMG_Epochs_R_resp =EMG_Epochs_R_resp[:,1,:]        
     
        EMG_Epochs_av_resp = np.concatenate((EMG_Epochs_L_resp,EMG_Epochs_R_resp)).mean(axis = 0)
        
        # Compute baseline --------------------------------------------
        base = np.mean(EMG_Epochs_av[:t0])
        EMG_Epochs_av = EMG_Epochs_av-base             #substract baseline from signal EMG locked
        EMG_Epochs_av_resp = EMG_Epochs_av_resp-base   #substract baseline from signal EMG response-locked
        
        ind_av_EMG_All[0, i,s,:] = EMG_Epochs_av    
        ind_av_EMG_All_resp[0,i,s,:] = EMG_Epochs_av_resp  



        # 2000ms stimulus duration 
        # EMG-locked -----------
        Epochs_L = MARK_seg_modified_continuous.segment(codeT0=locking_events_lag2000[i,0], tmin=min_stim, tmax=max_stim)
        EMG_Epochs_L = evt.getDataEpochs(EMG_data, Epochs_L)#(48, 2, 412)
        EMG_Epochs_L =EMG_Epochs_L[:,0,:]    
         
        Epochs_R = MARK_seg_modified_continuous.segment(codeT0=locking_events_lag2000[i,1], tmin=min_stim, tmax=max_stim)
        EMG_Epochs_R = evt.getDataEpochs(EMG_data, Epochs_R)    
        EMG_Epochs_R =EMG_Epochs_R[:,1,:]    

        EMG_Epochs_av = np.concatenate((EMG_Epochs_L,EMG_Epochs_R)).mean(axis = 0)        
    
        # response-locked -----------
        Epochs_L_resp = MARK_seg_modified_continuous.segment(codeT0=resp_locking_events_lag2000[i,0], tmin=min_resp, tmax=max_resp)
        EMG_Epochs_L_resp = evt.getDataEpochs(EMG_data, Epochs_L_resp)
        EMG_Epochs_L_resp =EMG_Epochs_L_resp[:,0,:]
        
        Epochs_R_resp = MARK_seg_modified_continuous.segment(codeT0=resp_locking_events_lag2000[i,1], tmin=min_resp, tmax=max_resp)
        EMG_Epochs_R_resp = evt.getDataEpochs(EMG_data, Epochs_R_resp)    
        EMG_Epochs_R_resp =EMG_Epochs_R_resp[:,1,:]        
     
        EMG_Epochs_av_resp = np.concatenate((EMG_Epochs_L_resp,EMG_Epochs_R_resp)).mean(axis = 0)
        
        # Compute baseline --------------------------------------------
        base = np.mean(EMG_Epochs_av[:t0])
        EMG_Epochs_av = EMG_Epochs_av-base            #substract baseline from signal EMG locked
        EMG_Epochs_av_resp = EMG_Epochs_av_resp-base  #substract baseline from signal EMG response-locked
        
        ind_av_EMG_All[1, i,s,:] = EMG_Epochs_av    
        ind_av_EMG_All_resp[1,i,s,:] = EMG_Epochs_av_resp  
        
    print("  >> Subject", s+1,"completed")



# Save data -----------------------------------------------------------------------------------
df_PRB_analysis_exp4 = df_PRB_analysis_exp4.drop(0) # delete 1sr row only for pandas creation
df_PRB_analysis_exp4.to_csv('df_PRB_analysis.csv', index=False) 
array_data = df_PRB_analysis_exp4.values

if saving:
    np.save('df_PRB_analysis_new.npy', array_data)
    np.save('ind_av_EMG_All.npy', ind_av_EMG_All)    
    np.save('ind_av_EMG_All_resp.npy', ind_av_EMG_All_resp)     
    np.save('data_export.npy', export)

print('\n>>Completed without error.')
# [X]