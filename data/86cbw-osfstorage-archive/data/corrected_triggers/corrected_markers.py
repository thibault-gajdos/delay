"""
Corrected markers
--> for technical reasons, RTs collected in the behavioral datafile are more accurate than the response triggers in the EMG file.
--> for this reason, let's adjust the response triggers based on the behavioral datafile (you can ignore this section, as it is a bit technical). 

"""

# Libraries 
import numpy as np               # To manipulate array table
import mne                       # MNE python library to process EEG (and also MEG) data. You can check MNE version with 'print(mne.__version__)'
import os                        # Useful for directory specification later on
import sys                       # To add a folder to the python path
import copy                      # To copy your python object
import matplotlib.pyplot as plt  # Matplotlib is a plotting library
import pandas as pd              # pandas is a dataframe library
import arviz as az               # arviz is a library for plotting 
current = os.path.join(os.getcwd()) 
os.chdir(r"C:\Users\erc\iCloudDrive\ACHARD_these\EMG_method\toolbox_EMG")
import emgTools as emg
import evtTools as evt
from Events import Events
import useMne
import vizEMG

# Plot parameters
az.style.use('arviz-darkgrid') # style fonts for plots
np.set_printoptions(5, suppress=True) # display 5-decimal numbers

# Path ------------------------------------------
path = os.path.join(current, "data_saved")
path_bdf = os.path.abspath(os.path.join(current, "..", "data", "bdf"))
path_behavior = os.path.abspath(os.path.join(current, "..", "data", "behavioral"))
path_corrected_triggers = os.path.abspath(os.path.join(current, "..", "data", "corrected_triggers"))
os.chdir(path)


# Select subject --------------------------------------
subject = 1

# Load bdf file ---------------------------------------
raw = mne.io.read_raw_bdf('s'+str(subject) + ".bdf", preload = True)
sf = raw.info['sfreq']
print ("successful EMG data loading")

# Extract trigger events from raw files
raw_events = mne.find_events(raw) 
#first column is the sample of the event (1024 samples = 1 second) since recording onset, third column is trigger (either stimulus or response)

# Load behavioral data --------------------------------
behavioral_data = np.load(path_behavior + f'\s{subject}.npy') 
RT = behavioral_data[:, 12].astype(np.float)*1000   #extract RTs and convert in milliseconds 
acc = behavioral_data[:, 11].astype(np.float)
direction = behavioral_data[:, 7].astype(np.float)


# If technical problem (delete after), e.g., trial no saving in bdf for s7 and s8, and s17 because 1 frame (16 ms) anticipation at RS (RT = 11ms).
if subject == 7:
    sample_trial = (raw_events[767,0] + (RT[197]/1000)*sf) +2
    raw_events = np.insert(raw_events, 768, [sample_trial, 0, 100], axis = 0)
if subject == 8:
    sample_trial = (raw_events[706,0] + (np.abs(RT[176])/1000)*sf)
    raw_events = np.insert(raw_events, 707, [sample_trial,0, 101], axis = 0)
if subject == 17:
    sample_trial = (raw_events[2596,0] + (np.abs(RT[650])/1000)*sf)
    raw_events = np.insert(raw_events, 2597, [sample_trial,0, 101], axis = 0) # left correct  (see on .npy)
if subject == 18:
    sample_trial = (raw_events[2875,0] + (np.abs(RT[720])/1000)*sf)
    raw_events = np.insert(raw_events, 2876, [sample_trial,0, 103], axis = 0) # right correct 


# First let's plot the difference between RTs from the behavioral datafile and the bdf file to see the difference ---
# extract RT fom the EMG/EEG (.bdf) file
RT_bdf = []
RT_without_misses = []#select RT from the behavioral file (removing missed trials)
behav_info = [] # conserve behvioral info
counter = 0# to keep track of each trial
premature_response = 0 #count number of premature response

vector_bad_experiment_code_trials = []
# 0 = non valid (premature response, missed trial, lag >= 20ms)  /  1= valid trial (lag < 20 ms)
miss = 0
check = []
for i in range (len(raw_events)):#loop over triggers from the EMG file
    if raw_events[i, 2] == 250: #250 is the trigger for 'no response' (missed trial)
        counter = counter + 1 
        vector_bad_experiment_code_trials.append(0)
        miss += 1
    
    if raw_events[i, 2] in [99, 229, 231]: # e.g., note technical error [7295434, 135, 231] for raw_events[3350] for s14
        counter = counter + 1 
        vector_bad_experiment_code_trials.append(0)
        print("errors:",i)

    if raw_events[i, 2] in [100, 101, 102, 103]:# a response has been given by the subject at right time
        # convert the sample corresponding to each event  (stimulus and response) into time (in ms); 
        # to do that, we divide the sample by the sampling rate (1024) and multiply the result by 1000
        # we then compute time of the response - time of the stimulus to get the RT
        a = (raw_events[i, 0]/sf)*1000 - (raw_events[i-1, 0]/sf)*1000  # RT bdf file
        b = RT[counter] # RT behavior data
        RT_bdf.append(a)
        RT_without_misses.append(b)
        behav_info.append(behavioral_data[counter,(1,2, 4, 11,12,13)])
        
        # Note if difference between RT behav. We consider that the trial is valid, otherwise, we consider that the trial is not valid (because of a lag > 20 ms between the two files, which is too much for a RT)
        if (a - b) < 20:
            vector_bad_experiment_code_trials.append(1)
            check.append(1)
        if (a - b) >= 20:
            vector_bad_experiment_code_trials.append(0)
            check.append(0)
            print(i, counter)  
        counter = counter + 1

    if raw_events[i, 2] in [200,201,202,203]: # a premature response
        vector_bad_experiment_code_trials.append(0)
        premature_response += 1
        counter += 1


# Convert list to numpy array 
RT_bdf = np.array(RT_bdf) 
RT_without_misses = np.array(RT_without_misses)
behav_info = np.array(behav_info)

# Compute difference between RT in the bdf file and RT in the behavioral file  -------------------
lag = RT_bdf - RT_without_misses 
check = np.array(check)

vector_bad_experiment_code_trials = np.array(vector_bad_experiment_code_trials)
print ("Bad trials with a lag > 20 ms, missed trial or premature response:",len(vector_bad_experiment_code_trials[vector_bad_experiment_code_trials==0]))
print("Premature response :", premature_response,"/960 =", round((premature_response/960*100),2), "%")


# Plot lag before processing -----------------------------------------------------------------------
fig, axs = plt.subplots(1, 3,   figsize=(15,5), dpi=300)
x = np.arange(1, len(lag)+1)
y = lag
c = behav_info[:, 1]
data = pd.DataFrame({"Trial": x, "Lag": y, "Category": c})
groups = data.groupby("Category")
plt.subplot(1,3,1)
for name, group in groups:
    plt.plot(group["Trial"], group["Lag"], marker="o", linestyle="", label=name)
plt.legend(bbox_to_anchor=(-0.09, 0.8),  borderaxespad=0, title = 'Block')
plt.title("Lag before without premature response")

# Lag after trimming
lag_trim = lag[check==1]
plt.subplot(1,3,2)
plt.scatter(np.arange(1, len(lag_trim)+1), lag_trim)
plt.title("Lag after trimming bad trials\n(" + str(len(vector_bad_experiment_code_trials[vector_bad_experiment_code_trials==0])) + " with " +
          str(premature_response) + " premature responses)")





# Now correct the markers in the corrected_events variable ---------------------------------
RT_converted = np.round((1024*RT)/1000)      # convert RT from the behavioral file into samples
corrected_events = copy.deepcopy(raw_events) # make a deep copy of the event file
counter = 0
lag = []
for i in range(len(corrected_events)):
    if corrected_events[i, 2] == 250:
        counter = counter + 1 #don't modify the "missed trial" marker
    
    # technical error for trials in s14, 17, 18
    if corrected_events[i, 2] in [229,231]:
        corrected_events[i, 0] = corrected_events[i-1, 0] + RT_converted[counter]
        counter = counter + 1 

    if corrected_events[i, 2] in [200,201,202,203]:
        counter = counter + 1 #don't modify the "premature response" marker
        
    if corrected_events[i, 2] in [100,101,102,103]:
        if vector_bad_experiment_code_trials[counter] == 1:
            corrected_events[i, 0] = corrected_events[i-1, 0] + RT_converted[counter]
            #adjust the timing of response markers to match the behavioral datafile 

        if vector_bad_experiment_code_trials[counter] == 0:
            # 99 marker will be > 20 ms lag
            corrected_events[i, 2] = 99 #change marker to "bad trial"
        counter = counter + 1



# Compute corrected RTs and compare with behavioral datafile -------------------------------------
RT_behav_clean = []
RT_bdf_corrected = []
counter = 0
for i in range (len(corrected_events)):
    if corrected_events[i, 2] in [250, 200,201,202,203, 99]: 
        counter = counter + 1
    if corrected_events[i, 2] in [100,101,102,103, 229, 231]:
        RT_bdf_corrected.append((corrected_events[i, 0]/1024)*1000 - (corrected_events[i-1, 0]/1024)*1000)
        RT_behav_clean.append(RT[counter])
        counter = counter + 1


RT_bdf_corrected = np.array(RT_bdf_corrected)  
RT_behav_clean = np.array(RT_behav_clean)  
lag_after = np.array(RT_behav_clean - RT_bdf_corrected)     

plt.subplot(1,3,3)
plt.scatter(np.arange(1, len(lag_after)+1), lag_after, color = 'green')
plt.title("Lag after correction") #should be centered on 0
plt.show()


# Save corrected marker data -----------------------------------------------------------
np.save(path_corrected_triggers + f'\s{subject}_corrected_markers.npy', corrected_events)

print(">>> Succefull run, subject", subject)
# [X]