""" Average accuracy, RT, PMET MET """

# Importation and Load data
import sys
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from matplotlib import style
import arviz as az
import scipy.stats as st
np.set_printoptions(5, suppress=True)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
from pathlib import Path
import matplotlib.pyplot as plt
import arviz as az
current = os.path.join(os.getcwd()) 

path_vector_TT = os.path.abspath(os.path.join(current, "..", "..","data", "vector_trial_type"))
path_data_export = os.path.abspath(os.path.join(current, "..", "..","data", "data_export"))
df = np.load(path_data_export + "\data_export.npy")

## ---------- Initialise data containers 
subjects_id = np.unique(df[:, 0])
nb_subjects = len(subjects_id) 
coh_stimLag = 3*2  # 3mc x 2stim lag

mean_RT = np.zeros((nb_subjects, coh_stimLag))
mean_PMT = np.zeros((nb_subjects, coh_stimLag))
mean_MT =  np.zeros((nb_subjects, coh_stimLag))
mean_RT_all = np.zeros((nb_subjects, coh_stimLag))
mean_PMT_all = np.zeros((nb_subjects, coh_stimLag))
mean_MT_all =  np.zeros((nb_subjects, coh_stimLag))
accuracy = np.zeros((nb_subjects, coh_stimLag))

miss_RS = np.zeros(nb_subjects)
anticip_RS = np.zeros(nb_subjects)
premature_response_RS = np.zeros(nb_subjects)
rejet_RS = np.zeros((nb_subjects))

# Coherence and stimLag
possible_coh_levels = np.array([0, 0.12, 0.5 ])
possible_stimLag = np.array([0.3,2.])

# stat trial type
stat_behav = np.zeros((nb_subjects, 6)) # 21 subjects, 5 trials type


s = 0
for sub_base in subjects_id:
    sub_base = int(sub_base)
    data_suj = df[(df[:, 0] == sub_base)] #extract data from subject
    trial = data_suj[:, 1]
    coherence = data_suj[:, 2]
    acc = data_suj[:, 3]
    RT = data_suj[:, 4]
    PMT = data_suj[:, 5]
    MT = data_suj[:, 6]
    TT_behav = data_suj[:, 7] # behavior + emg trial type
        # (0 Response before Response Signal (RS), DR blocks only)
        # (1 Missed Trial)
        # (2 Anticipation trials)
        # (3 Could not detect EMG onsets)
        # (4 technical problem, remove this trial) #will override all other trial type
        # (5 Valid trial, no problem cited above)
    stimLag = data_suj[:, 14]
    before_RS = data_suj[:, 15] 

    
    # Load EMG/EEG artifact vector data by subject
    TT = np.load(path_vector_TT + "\s" + str(sub_base) + "_vector_TT_EMG_EEG.npy")
        # (0 Response before Response Signal (RS))
        # (1 Missed Trial)
        # (2 Anticipation RT < 200 ms) 
        # (3 Could not detect EMG onsets)
        # (4 technical problem) 
        # (5 Valid trial, to include)
        # (6 if EEG artifacts in valid trials)
    
    # compute stat trial type
    tot_trial = len(TT)
    for i in range(6):
        stat_behav[s,i] = (len(TT_behav[(TT_behav == i)])/ tot_trial)*100
    
    counter = 0                    
    for d in range(len(possible_stimLag)):
        # by coherence
        for c in range(len(possible_coh_levels)):
            mask = (TT == 5) & (coherence == possible_coh_levels[c]) & (stimLag == possible_stimLag[d]) & (before_RS == 1)
                
            mean_RT[s,counter] = np.mean(RT[mask & (acc==1)])
            mean_PMT[s,counter] = np.mean(PMT[mask & (acc==1)])
            mean_MT[s,counter] = np.mean(MT[mask & (acc==1)])
            
            mean_RT_all[s,counter] = np.mean(RT[mask & (acc != -1)])
            mean_PMT_all[s,counter] = np.mean(PMT[mask & (acc != -1)])
            mean_MT_all[s,counter] = np.mean(MT[mask &  (acc != -1)])
            
            accuracy[s, counter] = np.mean(acc[mask & (acc != -1) ])
        
            counter += 1
            
    s += 1
    print("Subject", sub_base,"completed")


# --- stat TT ---
column_names = [
    "Response before Response Signal (RS)",
    "Missed Trial",
    "Anticipation RT < 200 ms",
    "Could not detect EMG onsets",
    "Technical problem",
    "Valid trial, to include"
]

df_stat_behav = pd.DataFrame(np.round(stat_behav, decimals=2), columns=column_names)

#GA stat accross subjects (in %)
GA_stat = np.round(np.mean(stat_behav, axis = 0),2)
min_stat = np.round(np.min(stat_behav, axis = 0),2)
max_stat = np.round(np.max(stat_behav, axis = 0),2)

df = {
    "GA_stat_TT (%)": GA_stat,
    "min": min_stat,
    "max": max_stat
}
df_GA_stat = pd.DataFrame(df, index=column_names)




# Export data JASP ---------------------------------------------------------
id_subject = list(range(1, 23, 1))*6

stimLag = list([[300]*22*3,[2000]*2*3])
stimLag = [x for y in stimLag for x in y]

coherence = list([[0]*2,[.12]*2, [.50]*2]*3)
coherence = [x for y in coherence for x in y]


# Create data csv (acc = 1)
mean_RT = pd.DataFrame(mean_RT)
mean_RT.columns =  ['mean_RT_coh0_300','mean_RT_coh12_300','mean_RT_coh50_300',
                    'mean_RT_coh0_2000','mean_RT_coh12_2000','mean_RT_coh50_2000']


mean_PMT = pd.DataFrame(mean_PMT)
mean_PMT.columns = ['mean_PMT_coh0_300','mean_PMT_coh12_300','mean_PMT_coh50_300',
                    'mean_PMT_coh0_2000','mean_PMT_coh12_2000','mean_PMT_coh50_2000']


mean_MT = pd.DataFrame(mean_MT)
mean_MT.columns = ['mean_MT_coh0_300','mean_MT_coh12_300','mean_MT_coh50_300',
                    'mean_MT_coh0_2000','mean_MT_coh12_2000','mean_MT_coh50_2000']

mean_accuracy = pd.DataFrame(accuracy)
mean_accuracy.columns = ['accuracy_coh0_300','accuracy_coh12_300','accuracy_coh50_300',
                        'accuracy_coh0_2000','accuracy_coh12_2000','accuracy_coh50_2000']


# Create data csv (correct + incorrect - acc != -1)
mean_RT_all = pd.DataFrame(mean_RT_all)
mean_RT_all.columns =  ['mean_RT_coh0_300','mean_RT_coh12_300','mean_RT_coh50_300',
                    'mean_RT_coh0_2000','mean_RT_coh12_2000','mean_RT_coh50_2000']


mean_PMT_all = pd.DataFrame(mean_PMT_all)
mean_PMT_all.columns = ['mean_PMT_coh0_300','mean_PMT_coh12_300','mean_PMT_coh50_300',
                    'mean_PMT_coh0_2000','mean_PMT_coh12_2000','mean_PMT_coh50_2000']


mean_MT_all = pd.DataFrame(mean_MT_all)
mean_MT_all.columns = ['mean_MT_coh0_300','mean_MT_coh12_300','mean_MT_coh50_300',
                    'mean_MT_coh0_2000','mean_MT_coh12_2000','mean_MT_coh50_2000']




# Saving -----------------------------------
# mean_accuracy.to_csv("accuracy.csv",index=False)
# mean_RT_all.to_csv("mean_RT_all.csv",index=False)
# mean_PMT_all.to_csv("mean_PMT_all.csv",index=False)
# mean_MT_all.to_csv("mean_MT_all.csv",index=False)


print('Successing process without errors.')
# [X]