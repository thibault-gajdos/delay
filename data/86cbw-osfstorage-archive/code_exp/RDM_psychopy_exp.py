# -*- coding: utf-8 -*-
 
"""
Experiment: EMG EEG Experiment in delayed response 3s RDM task
@jachard

---- Columns output data ----

|  | Columns Name           | Signification      | Var   | Details                                                                                                        |
|  |                        |                    |       |                                                                                                                |
|0 | subject_id             | ID Subject         | Int   | Subject number.                                                                                                |
|1 | trial_tot              | Trial number tot   | Int   | Number of trial in total for all block (1 -> 960 trials) 150 trials/conditions                                 |
|2 | block                  | Block number       | Int   | Number of trial in total for all block.                                                                        |
|3 | block_type             | Type of block      | Str   | There is Free Response (FR) or Response Signal (RS - Delay response with 3sec delay).                          |
|4 | trial_block            | Trial number block | Int   | Number of trial for block (1 -> 96 trials)                                                                     |
|5 | condition_coherence    | Motion coherence   | Int   | Motion coherence of trial. Condition of coherence possible: .0 .12  0.5                                        |
|6 | stim_pres              | Time of pres. stim | Float | Time of presentation of stimulus during stimulus period: 300 ms or 2000ms, -1 if FR                            |
|7 | direction_motion       | Dots direction     | Int   | Direction of dots in RDM for the trial. Left is 180, Right is 0.                                               |
|8 | key_required           | Expected key       | Int   | The key response expected for a correct response (1-left or 2-right).                                          |
|9 | key_pressed            | Pressedey          | Int   | The key response pressed (1-left or 2-right).                                                                  |
|10| response_moment        | Time of response   | Int   | The moment when the participant pressed a key can be too early (0), at the right time (1), or to late (-1).    |
|11| accuracy               | Accuracy           | Int   | The accuracy pf each subject for the task.                                                                     |
|12| rt                     | RT                 | Float | Response Time of participant from start of stimulus to button press (in ms), before the disapear of cue  .     |
|13| trigg_stim             | Stimulus trigger   | Int   | The stimulus trigger indicate conditions select of trial (delay, coherence, direction of dots). See * below.   |
|14| trigg_resp             | Response trigger   | Int   | The response trigger indicate response subject (correct or not, left or right). See ** below.                  |




---- TRIGGERS -----
- Fixing Cross  : 2
- Stimulus Onset : see Table below stimulus period
- Stimulus Offset (delay period- 3s) : 1
- Response Signal (response period) : trigger stimulus Onset + 100
- Response : 100 - 101 - 102 - 103 if response after response signal, or 200 - 201 - 202 - 203 if before response signal, or 250 if too late 


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


"""


from psychopy import core, data, event, visual, gui, monitors 
from psychopy.tools.monitorunittools import deg2pix
import numpy as np
from scipy.ndimage import maximum_filter
from scipy import sparse ## for compressing matrices
import os
from psychopy.hardware import keyboard
from psychopy.preferences import prefs
prefs.hardware['audioLib'] = ['PTB']
from psychopy import sound
import serial
import time 


#-----------------------------------  Variables ----------------------------------------
###### Activate port for toolkits for EMG
activate_port = False
port = serial.Serial("COM5", baudrate=115200) if activate_port else None


#-------------------------  RDM Parameters  ------------------------------- 
expName = "EEMG_delay_"
pc_manip = True

# >>>>>> Change here Folder <<<<<<< 
if pc_manip:
    folder = r"C:\Users\mservant\Desktop\Johan\RDM_EEG_EMG_delay\data_10_50\csv"
    image_path = r"C:\Users\mservant\Desktop\Johan\RDM_EEG_EMG_delay\oeil.jpg"
    folder_dots = r"C:\Users\mservant\Desktop\Johan\RDM_EEG_EMG_delay\data_10_50\csv\folder_dots"
else:
    folder = r"C:\Users\erc\iCloudDrive\ACHARD_these\EEMG_delay\experiment\test\csv"
    image_path = r"C:\Users\erc\iCloudDrive\ACHARD_these\EEMG_delay\experiment\oeil.jpg"



# Choose VDs
coherence_levels = np.array([0., .12, .50]) # 3 motions coherences levels 
stim_pres_level = np.array([18, 120]) # 2 duration of stimulus presentation 300ms and 2000ms in frame



# Choose globals variables experiment
touchIs=['1','2']    # response touchs
deadline = 4.        # max time to make a response 4sec
delay = 180          # delayed response fixe to 3 seconds in frame

# Screen Frame
jitter = np.array([36, 48, 60]) # 600ms, 800ms, 1000ms
eye_time = 90 # Pause ITI for flashing eyes 1.5s

repeat_stim_type_blockRS = 8 # Put '8' - Enter 8 for 96 trials/block, 12 possibles stim (2 directions left/right x 3 coherences x 2 stim duration 0.3/2s) - Participant see 16 times stimulus type by block (16 time, 6 possible stimulus )
number_block =  10           # Put '10' for 10 blocks of 96 trials (total 960 trials) - 160 trials/conditions


# Do you want Train ?
train = True          # To do train or not
nb_repeat_train = 3    # Enter 3 for 12 possible stim, repeat 3 times random  (36 trainning trials) 18 with feedback, 18 without feedback
 



#-------------------------- Experiment infos -------------------------------
# Experiment infos data 
expInfo = {'participant': ''}
date_exp = data.getDateStr() 
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, order=['participant'])
if not dlg.OK: core.quit() #user pressed cancel

filename = expName + '_s' + expInfo['participant'] + '_' + date_exp
running_filename = folder + '/' + filename + '_TEMP.txt'
output = "subject_id, trial_tot, block, block_type, trial_block, condition_coherence, stim_pres, direction_motion, key_required, key_pressed, response_moment, accuracy, RT, trigg_stim, trigg_resp\n"
with open(running_filename, 'a') as f_handle: f_handle.write(output)  # with open close automatically

subject_id = int(expInfo['participant'])

#---------------------------- Window parameter ---------------------------------
if pc_manip:
    win = visual.Window(fullscr=True, monitor='mathieu', size=(1920, 1080), color=[0,0,0], colorSpace='rgb255', units='deg')
else:
    win = visual.Window(fullscr=True, monitor='mathieu', size=(1920, 1080), color=[0,0,0], colorSpace='rgb255', units='deg')
police = 'Consolas' # monospaced font
win.setMouseVisible(False) # hide the mouse cursor!!!

# font size and settings
size_instruc = .5
size_feedback = .7

text_dict = {'win':win, 'ori':0, 'font':police, 'pos':[0, 0], 'colorSpace':'rgb255'}
instruc_dict = {**text_dict, 'height':size_instruc, 'color':[255,255,255]}
feedback_dict = {**text_dict, 'height':size_feedback} 
cross_dict  = {**text_dict, 'height':1.3} 
dot_dict = {**text_dict, 'height':5} 



#-----------------------------------------------------------------------------------------------
#  ________   __    __   __    __    ________   ________   __   ________   __    __   ________ 
# |   _____| |  |  |  | |   \ |  |  |   _____| |__    __| |  | |   __   | |   \ |  | |   _____|
# |  |___    |  |  |  | |    \|  |  |  |          |  |    |  | |  |  |  | |    \|  | |  |_____
# |   ___|   |  |  |  | |  |\    |  |  |          |  |    |  | |  |  |  | |  |\    | |_____   |
# |  |       |  |__|  | |  | \   |  |  |_____     |  |    |  | |  |__|  | |  | \   |  _____|  |
# |__|       |________| |__|  \__|  |________|    |__|    |__| |________| |__|  \__| |________|

#-----------------------------------  Programs functions ----------------------------------------
#------------------------------------------------------------------------------------------------
def instruction_exp(instruction, touchIs):
    kb_delay = 0
    win.callOnFlip(resp.clearEvents, eventType='keyboard') 
    while True:  
        if kb_delay:#check keyboard buffer AFTER first draw of stimulus and clock reset
            theseKeys = resp.getKeys(keyList=touchIs, waitRelease=False) 
            if len(theseKeys):
                break
            if defaultKeyboard.getKeys(keyList=["escape"], waitRelease=False):
                core.quit()   
        instruction.draw()
        win.flip()
        if not kb_delay:
            kb_delay = 1
    win.flip()
    core.wait(1) 
    return None


def find_response(trigg_stim, touchIs):
    global trigg_right, trigg_left
    # Find correct response
    if trigg_stim in trigg_left:
        correct_resp = touchIs[0]
        direction_trial = 180 #180 is leftward
    elif trigg_stim in trigg_right:
        correct_resp = touchIs[1]
        direction_trial = 0  #0 is rightward
    else: raise ValueError('Stimulus trigger not in the list...')

    # Set experiment condition
    # Coherence
    if trigg_stim in trigg_coh0: coherence_trial = coherence_levels[0]
    elif trigg_stim in trigg_coh10: coherence_trial = coherence_levels[1]
    elif trigg_stim in trigg_coh40: coherence_trial = coherence_levels[2]
    else: raise ValueError('No motion coherence was found in response...')

    # Stim presentation time
    if trigg_stim in trigg_RS_300: stim_pres_trial = stim_pres_level[0]
    elif trigg_stim in trigg_RS_2000: stim_pres_trial = stim_pres_level[1]
    else: raise ValueError('No duration of stimulus was found in response...')

    return correct_resp, direction_trial, coherence_trial, stim_pres_trial

"""
|                 | Stim Pres |  coherence  0.  |  coherence .10  |  coherence .40  |  
|                 |           |                 |                 |                 |  
| RS 300          |  300 ms   |   L 21 / R 24   |   L 22 / R 25   |   L 23 / R 26   |  
| RS 2000         |  2000 ms  |   L 31 / R 34   |   L 32 / R 35   |   L 33 / R 36   |  
"""


#------------------------------------------------------------------------------------------------------------------------------------------------------------
#  __   __    __   ________   ________   ________    __    __   ________   ________   __   ________   __    __   ________ 
# |  | |  \  |  | |   _____| |__    __| |   __   |  |  |  |  | |   _____| |__    __| |  | |   __   | |   \ |  | |   _____|
# |  | |    \|  | |  |_____     |  |    |  |__| _|  |  |  |  | |  |          |  |    |  | |  |  |  | |    \|  | |  |_____
# |  | |  |\    | |_____   |    |  |    |   __  \   |  |  |  | |  |          |  |    |  | |  |  |  | |  |\    | |_____   |
# |  | |  | \   |  _____|  |    |  |    |  |  \  \  |  |__|  | |  |_____     |  |    |  | |  |__|  | |  | \   |  _____|  |
# |__| |__|  \__| |________|    |__|    |__|   \__\ |________| |________|    |__|    |__| |________| |__|  \__| |________|

#-----------------------------------------------------  Instructions text ------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------


instructions_position = visual.TextStim(**instruc_dict,
        text =u"Bienvenue ! Prenez la commande gauche (celle associée à un sparadrap blanc) avec votre main gauche, et la commande droite avec votre main droite.\
 Posez vos mains sur le coussin, puis placez vos pouces juste au dessus de chaque bouton.\
Ne croisez jamais les mains. Ne prenez jamais la manette gauche avec la main droite, et la manette droite avec la main gauche.\
Détendez-vous, et relâchez le plus possible les muscles de vos mains tout en gardant cette position.\
Vous devrez essayer de conserver cette position et ce relâchement pendant toute la durée de l'expérience (sauf pendant les pauses), car les contractions musculaires (en dehors des réponses) pertubent le signal EMG. \n\n\
Par ailleurs, vous devrez toujours regarder la croix de fixation au centre de l'écran. \
Evitez le plus possible de bouger et cligner des yeux, car les mouvements oculaires perturbent le signal électrique cérébral (EEG). \
Vous pourrez uniquement bouger et cligner des yeux entre chaque essai (si vous le souhaitez), lorsque le dessin d'un oeil apparaîtra. \n\n\
Prenez le temps de bien relire ces instructions, car elles sont très importante pour la réussite de l'expérience. \n\n\
Pressez n'importe quel bouton pour découvrir les instructions de l'expérience.") 


instructions_RDM_task = visual.TextStim(**instruc_dict,
        text =u"Dans chaque essai, un ensemble de points blancs en mouvement sera présenté sur un fond noir. \
Une partie des points se déplacera vers la gauche ou vers la droite, et l'autre partie se déplacera aléatoirement. \
Vous devrez presser le bouton gauche avec votre pouce gauche si les points se dirigent vers la gauche, pressez le bouton droit avec votre pouce droit si les points se dirigent vers la droite.\n\n\
Si de nombreux points se dirigent vers la gauche ou vers la droite, l'essai sera relativement facile. Si très peu de points se dirigent vers la droite ou vers la gauche, l'essai sera très difficile. \
Par conséquent, il est tout à fait normal d'avoir une forte incertitude sur la réponse à fournir dans certains essais. Nous vous demandons toutefois de fournir une réponse à chaque essai, même si vous avez parfois l'impression de répondre au hasard.\n\n\
Pressez sur n'importe quel bouton pour continuer.")  


instructions_RS_training_0 = visual.TextStim(**instruc_dict,
        text =u"Dans cette expérience, vous devrez dire si les points en mouvement vont à droite ou à gauche en appuyant sur les boutons de réponses associés, après un temps d'attente imposé.\n\
Concretement pour chaque essai vous verrez une croix de fixation bleue suivi des points en mouvement, et à nouveau la croix bleue. Après un certain temps la croix deviendra jaune.\
Vous devrez donner votre réponse uniquement lorsque que la croix bleue devient jaune.\n\n\
Fixer toujours la croix central et cligner des yeux entre chaque essai lorsque le dessin d'un oeil apparaîtra.\n\n\
Vous allez d'abord effectuer quelques essais d'entraînement pour vous familiariser avec cette tâche.\
Dans ces premiers essais, nous vous indiquerons si vous avez fourni ou non la bonne réponse.\n\n\
Pressez sur n'importe quel bouton pour débuter l'entrainement.")  

instructions_RS_training_1 = visual.TextStim(**instruc_dict,
        text =u"Félicitation ! Vous avez découvert l'entrainement.\n\n\
Vous allez continuer l'entrainement sur cette tâche, mais maintenant vous n'aurez plus de feedback sur vos réponses.\n\n\
Pressez sur n'importe quel bouton pour continuer.")  


instructions_RS_training_2 = visual.TextStim(**instruc_dict,
        text =u"Félicitation ! Vous avez réussi l'entrainement.\n\n\
Rappelez-vous que dans cette tâche, vous devrez répondre le plus rapidement et précisement possible uniquement lorsque la croix de fixation devient jaune après le temps d'attente.\n\n\
Pressez sur n'importe quel bouton pour commencer l'expérience.")  



instructions_testing = visual.TextStim(**instruc_dict,
        text =u"Vous allez maintenant effectuer 10 blocs de 96 essais. Chaque bloc d'essais sera séparé par une petite pause.\n\
Vous ne recevrez pas de feedback correct/incorrect après vos réponses.\n\n\
Fixer toujours la croix central et cligner des yeux uniquement entre chaque essai devant l'écran correspondant.\n\n\
Merci de prévenir l'expérimentateur que vous allez débuter l'expérience.\n\n\
Appuyez sur n'importe quel bouton pour commencer.")   


instructions_pause = visual.TextStim(**instruc_dict,
        text =u"Prenez une petite pause si vous en ressentez le besoin, puis pressez n'importe quel bouton pour afficher le type de bloc suivant.\n\n\
IMPORTANT : si vous avez lâché les commandes, prenez soin de bien reprendre la commande gauche (celle avec le sparadrap blanc) avec votre main gauche, et la commande droite avec votre main droite.")  


fin=visual.TextStim(**instruc_dict,
    text=u"Fin de l'expérience. Merci beaucoup pour votre participation.")


# Feedback for training
correct_decision = visual.TextStim(**feedback_dict,
        text =u"CORRECT", color=[255,255,255])
        
incorrect_decision = visual.TextStim(**feedback_dict, 
        text =u"INCORRECT", color=[255,255,255])

trop_tard =visual.TextStim(**feedback_dict,
        text =u"TROP TARD !\n\nRépondez plus vite s'il vous plaît.", color=[255,255,255])

trop_tot = visual.TextStim(**feedback_dict,
        text =u"TROP TOT !\n\nRépondre lorsque la croix centrale devient jaune !", color=[255,255,255])



#------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------- Initialisation RDM --------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
### Initialize keyboards
resp = keyboard.Keyboard()
defaultKeyboard = keyboard.Keyboard()#for escaping only

####################### Parameters for RDK stimuli #######################
framerate = 60.0
speed_deg_sec = 10.5 #in degrees/s; 12 in Pilly in Seitz
dot_density = 16.7 #in deg-2 s-1
rayon_cercle = 9.0 #in deg - in experimental computer we mesured 15.42cm for diameter correspond to 11.78° apertune circular at 75cm from screen
number_of_dots = int(np.round(dot_density * np.pi * np.square(rayon_cercle) * (1/framerate)))
dotsize = 4

##########################################################################

### Dots creation 
dots = visual.DotStim(win=win, name='random_dots',
    nDots= number_of_dots, 
    dotSize=dotsize,#in pixels
    units = 'deg',
    speed= speed_deg_sec/framerate, #in degrees per frame
    dir= 0, #in degrees  --> manipulated
    coherence= 0.0,    # --> manipulated
    fieldPos=[0.0, 0.0], 
    fieldSize=rayon_cercle*2, 
    fieldShape='circle',
    signalDots='different', # on each frame the dots constituting the signal could be the same as on the previous frame or different. If 'same', participants may determine the direction of dots based on a single dot.
    noiseDots='position',   #can be set to 'direction' or 'location'. 'location' has the disadvantage that the noise dots not only have a random direction but alos a random speed (whereas signal dots have a constant speed and constant direction)
    dotLife= -1,            #number of frames each dot lives. 
    color=[255,255,255], colorSpace='rgb255') #color of the dots

# Fixation Cross
fixation_cross_blue = visual.TextStim(**cross_dict, text=u"+", color=[0, 229, 255])
fixation_cross_yellow = visual.TextStim(**cross_dict, text=u"+", color=[255, 255, 98])



# Load image of eye
eye_image = visual.ImageStim(win, image=image_path, pos=(0, 0))
eye_image.size = (8, 5) # Reduce size




#------------------------ Experiment Conditions ----------------------------
### Response condition for each hand
trigg_RS = [21, 22, 23, 24, 25, 26,   31, 32, 33, 34, 35, 36]
trigg_RS_300 = [21,22,23,24,25,26]
trigg_RS_2000 = [31,32,33,34,35,36]

trigg_left =  [21,22,23, 31,32,33]
trigg_right = [24,25,26, 34,35,36]

trigg_coh0 =  [21,24, 31,34]
trigg_coh10 =  [22,25, 32,35]
trigg_coh40 =  [23,26, 33,36]


jitter_trigg = 2
stim_offset_trigg = 1

"""
|    Coherence    | Stim Pres |  coherence  0.  |  coherence .10  |  coherence .40  |  
|                 |           |                 |                 |                 |  
| RS 300          |  300 ms   |   L 21 / R 24   |   L 22 / R 25   |   L 23 / R 26   |  
| RS 2000         |  2000 ms  |   L 31 / R 34   |   L 32 / R 35   |   L 33 / R 36   |  
"""

# Insruction
instruction_exp(instructions_position, touchIs)
instruction_exp(instructions_RDM_task, touchIs)



#------------------------------------------------------------------------------------------------------------------------------------------------------------
#  ________   ________    ________   __   __    __   __   __    __   ________
# |__    __| |   __   |  |   __   | |  | |   \ |  | |  | |   \ |  | |   _____|
#    |  |    |  |__| _|  |  |__|  | |  | |    \|  | |  | |    \|  | |  | ____
#    |  |    |   __  \   |   __   | |  | |  |\    | |  | |  |\    | |  ||_   |
#    |  |    |  |  \  \  |  |  |  | |  | |  | \   | |  | |  | \   | |  |__|  |
#    |__|    |__|   \__\ |__|  |__| |__| |__|  \__| |__| |__|  \__| |________|

#-----------------------------  Training sequence -----------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
if train:
    
    instruction_exp(instructions_RS_training_0, touchIs)

    # Create training sequence
    block_RS = np.tile(trigg_RS, nb_repeat_train) # if nb_repeat_train = 3, 12 x 3 = 36 training trials 
    np.random.shuffle(block_RS)
    training_seq = np.vstack((block_RS[0:18],block_RS[18:36]))

    feedback = True

    for block in training_seq:
        trial_counter = 1
        for trigg_stim in block:

            # Find each variable necessary for trial
            correct_resp, direction_trial, coherence_trial, stim_pres_trial = find_response(trigg_stim, touchIs)
            
            # Update dots class
            dots.setFieldCoherence(coherence_trial)
            dots.setDir(direction_trial) 

            # Initialize variable
            key_pressed = -1
            rt = -1   # RT
            acc = -1  # accuracy
            response_moment = -1  # period of response, 0 if too early (before indice cue), 1 if at time (between indice cue and deadline, or FR), -1 if no response or > deadline

            # Initiate data contains dots coordinates
            dots_coord = []
            win.nDroppedFrames = 0 # to check frame
            win.recordFrameIntervals = True

            # warm up GPU
            for _ in range(3): win.flip()
            
            # --------- Display the eye picture for 1.5 seconds ---
            for _ in range(eye_time):
                eye_image.draw()
                win.flip()  

            # ------------ Display cross fixation (jitter time, random previously)
            if activate_port is True: port.write(bytes(jitter_trigg)) # trigger fixation cross appear in baseline screen 
            time_jitter = np.random.choice(jitter, 1)[0]
            for _ in range(time_jitter):
                fixation_cross_blue.draw()
                win.flip()  


            # Reset Clock and send trigg_stim to second pc at next win.flip()
            win.callOnFlip(resp.clock.reset)
            if activate_port is True: win.callOnFlip(port.write, bytes([int(trigg_stim)]))
            kb_delay = 0 

            ### Start Trial
            if (trigg_stim in trigg_RS):
                Resp_Signal_ON = 0   # 0 if reponse signal OFF (we are in stimulus or delay period), 1 if response signal ON (disapear of circle after stimulus and delay, we are in response period)
                Stim_ON = 1          # 0 if stimulus off, 1 if on
                
                # Time dipslay of experiment
                time_stim_delay = stim_pres_trial + delay # in frames

                count_frame = 0

                # Stimulus ONSET, delay and response display 
                while True:
                    if kb_delay: #check keyboard buffer AFTER first draw of stimulus and clock reset
                        
                        # Time clock (start with stimulus onset)
                        t = resp.clock.getTime()

                        # Stop trial deadline
                        if (Resp_Signal_ON==1) and (t > deadline):
                            trigg_resp = 250
                            if activate_port is True: port.write(bytes([trigg_resp]))
                            trop_tard.draw()
                            win.flip()
                            core.wait(2)
                            break

                        if (Resp_Signal_ON==1) and (t <= deadline):
                            fixation_cross_yellow.draw() # RS

                        # Touch response
                        theseKeys = resp.getKeys(keyList=touchIs, waitRelease=False) # response autorise
                        if len(theseKeys):
                            theseKeys = theseKeys[0] 
                            key_pressed = theseKeys.name

                            # Early response
                            if (Resp_Signal_ON == 0):
                                response_moment = 0
                                if key_pressed == correct_resp:
                                    acc = 1
                                    if key_pressed == touchIs[0]: trigg_resp = 201 # correct left button response
                                    if key_pressed == touchIs[1]: trigg_resp = 203 # correct right button response 

                                # Incorrect response
                                elif key_pressed != correct_resp:
                                    acc = 0
                                    if key_pressed == touchIs[0]: trigg_resp = 200 #incorrect left button response
                                    if key_pressed == touchIs[1]: trigg_resp = 202 #incorrect right button response 

                                rt = theseKeys.rt
                                if activate_port is True: port.write(bytes([trigg_resp]))                            

                                trop_tot.draw()
                                win.flip()
                                core.wait(2) 
                                break
                            
                            # Right time (after RS)
                            if (Resp_Signal_ON==1): 
                                
                                response_moment = 1
                                
                                if key_pressed == correct_resp:
                                    acc = 1
                                    if key_pressed == touchIs[0]: trigg_resp = 101 # correct left button response
                                    if key_pressed == touchIs[1]: trigg_resp = 103 # correct right button response 

                                # Incorrect response
                                elif key_pressed != correct_resp:
                                    acc = 0
                                    if key_pressed == touchIs[0]: trigg_resp = 100 #incorrect left button response
                                    if key_pressed == touchIs[1]: trigg_resp = 102 #incorrect right button response 

                                rt = theseKeys.rt 
                                if activate_port is True: port.write(bytes([trigg_resp]))

                                if feedback:
                                    feedback = correct_decision if key_pressed == correct_resp else incorrect_decision    
                                    win.flip()     
                                    feedback.draw()
                                    win.flip()# ? 
                                    core.wait(1.)  

                                win.flip()
                                break                        
                            
                        # 'Escp' to quit
                        if defaultKeyboard.getKeys(keyList=["escape"], waitRelease=False): core.quit()

                        # Draw circle and dots
                        if (Resp_Signal_ON == 0):
                            if (count_frame <= time_stim_delay):
                                fixation_cross_blue.draw() # draw circle during stim and delay period
                                if (count_frame <= stim_pres_trial): 
                                    dots.draw()
                                    #core.wait(1)

                        
                        # Trigger end of presentation dots
                        if (Resp_Signal_ON == 0) and (Stim_ON == 1) and (count_frame > stim_pres_trial): 
                            if activate_port is True: win.callOnFlip(port.write, bytes([int(1)]))
                            Stim_ON = 0
                    
                        if (Resp_Signal_ON == 0) and (count_frame > time_stim_delay) :
                            fixation_cross_yellow.draw() # RS
                            trigg_stim_delay = int(trigg_stim) + 100

                            win.callOnFlip(resp.clock.reset) # reset clock at the end of delay
                            if activate_port is True: win.callOnFlip(port.write, bytes([trigg_stim_delay]))       
                            Resp_Signal_ON = 1

                    if not kb_delay:
                        # Draw first dots
                        dots.draw()
                        fixation_cross_blue.draw()

                        resp.clearEvents(eventType='keyboard')
                        kb_delay = 1 

                    count_frame += 1
                    win.flip()  # function resp.clock.reset synchronized with the frame flip with callOnFlip

        if feedback == True:
            instruction_exp(instructions_RS_training_1, touchIs)
        feedback = False



instruction_exp(instructions_RS_training_2, touchIs)



#------------------------------------------------------------------------------------------------------------------------------------------------------------
#  ________   ________   ________   ________   __   __    __   ________
# |__    __| |   _____| |   _____| |__    __| |  | |   \ |  | |   _____|
#    |  |    |  |__     |  |_____     |  |    |  | |    \|  | |  | ____
#    |  |    |   __|    |____    |    |  |    |  | |  |\    | |  ||_   |
#    |  |    |  |_____   _____|  |    |  |    |  | |  | \   | |  |__|  |
#    |__|    |________| |________|    |__|    |__| |__|  \__| |________|

#-----------------------------  Experiment sequence --------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
if activate_port is True: port.write(bytes([5])) # Biosemi PauseOff = 5, to saving auto


### ------------------ Create each block type and randomize it ------------------
# Generate random sequence
block_RS = np.tile(trigg_RS, repeat_stim_type_blockRS) # if repeat_stim_type_blockRS = 5, 90 trials/block
np.random.shuffle(block_RS)
general_seq =  block_RS
for n in range(number_block-1):
    np.random.shuffle(block_RS)
    general_seq = np.vstack((general_seq,block_RS))

# checking general sequence
# assert general_seq.shape[1] == 96, print('Error in creation of general sequence, no 96 trials per block')
# assert general_seq.shape[0] == 12, print('Error in creation of general sequence, no 12 blocks')





# ------------------ Draw instructions start of experiment ------------------
instruction_exp(instructions_testing, touchIs)

#initialize Thibault's section for recording dot coordinates and dropped frames
mon = monitors.Monitor('mathieu')
ppd = int(deg2pix(1, mon, correctFlat=False)) ## ppd
n_pix =  int(2 * rayon_cercle * ppd) #image size in pixels
grid = np.linspace(-rayon_cercle, rayon_cercle, n_pix) ##grid for mesh
xx, yy = np.meshgrid(grid, grid) ##mesh frame
ntrials = general_seq.shape[1]*number_block
dropped = np.empty(ntrials)
dropped.fill(-1)


block_counter = 1
trial_tot = 1
data_container = []
trials_frame = []

# Loops in blocks
for block in general_seq:

    trial_counter = 1
    for trigg_stim in block:
        # Find RS block, and name block for output
        if trigg_stim in trigg_RS_300: type_block = "300"
        elif trigg_stim in trigg_RS_2000: type_block = "2000"
        else: raise ValueError('Oups, trigg_stim not in block...')

        # Find each variable necessary for trial
        correct_resp, direction_trial, coherence_trial, stim_pres_trial = find_response(trigg_stim, touchIs)
        
        # Update dots class
        dots.setFieldCoherence(coherence_trial)
        dots.setDir(direction_trial) 

        # Initialize variable
        key_pressed = -1
        rt = -1   # RT
        acc = -1  # accuracy
        response_moment = -1  # period of response, 0 if too early (before indice cue), 1 if at time (between indice cue and deadline, or FR), -1 if no response or > deadline

        # Initiate data contains dots coordinates
        dots_coord = []
        win.nDroppedFrames = 0 # to check frame
        win.recordFrameIntervals = True

        # warm up GPU
        for _ in range(3): win.flip()
        
        # --------- Display the eye picture for 1.5 seconds ---
        for _ in range(eye_time):
            eye_image.draw()
            win.flip()  

        # ------------ Display cross fixation (jitter time, random previously)
        if activate_port is True: port.write(bytes(jitter_trigg)) # trigger fixation cross appear in baseline screen 
        time_jitter = np.random.choice(jitter, 1)[0]
        for _ in range(time_jitter):
            fixation_cross_blue.draw()
            win.flip()  


        # Reset Clock and send trigg_stim to second pc at next win.flip()
        win.callOnFlip(resp.clock.reset)
        if activate_port is True: win.callOnFlip(port.write, bytes([int(trigg_stim)]))
        kb_delay = 0 

        ### Start Trial
        if (trigg_stim in trigg_RS):
            Resp_Signal_ON = 0   # 0 if reponse signal OFF (we are in stimulus or delay period), 1 if response signal ON (disapear of circle after stimulus and delay, we are in response period)
            Stim_ON = 1          # 0 if stimulus off, 1 if on
            
            # Time dipslay of experiment
            time_stim_delay = stim_pres_trial + delay # in frames

            count_frame = 0  # count number of frame to be sure that correspond to time duration

            # Stimulus ONSET, delay and response display 
            while True:
                if kb_delay: #check keyboard buffer AFTER first draw of stimulus and clock reset
                    
                    # Time clock (start with stimulus onset)
                    t = resp.clock.getTime()

                    # Stop trial deadline
                    if (Resp_Signal_ON==1) and (t > deadline):
                        trigg_resp = 250
                        if activate_port is True: port.write(bytes([trigg_resp]))
                        trop_tard.draw()
                        win.flip()
                        core.wait(2)
                        break

                    if (Resp_Signal_ON==1) and (t <= deadline):
                        fixation_cross_yellow.draw() # RS

                    # Touch response
                    theseKeys = resp.getKeys(keyList=touchIs, waitRelease=False) # response autorise
                    if len(theseKeys):
                        theseKeys = theseKeys[0] 
                        key_pressed = theseKeys.name

                        # Early response
                        if (Resp_Signal_ON == 0):
                            response_moment = 0
                            if key_pressed == correct_resp:
                                acc = 1
                                if key_pressed == touchIs[0]: trigg_resp = 201 # correct left button response
                                if key_pressed == touchIs[1]: trigg_resp = 203 # correct right button response 

                            # Incorrect response
                            elif key_pressed != correct_resp:
                                acc = 0
                                if key_pressed == touchIs[0]: trigg_resp = 200 #incorrect left button response
                                if key_pressed == touchIs[1]: trigg_resp = 202 #incorrect right button response 

                            rt = theseKeys.rt
                            if activate_port is True: port.write(bytes([trigg_resp]))                            

                            trop_tot.draw()
                            win.flip()
                            core.wait(2) 
                            break
                        
                        # Right time (after RS)
                        if (Resp_Signal_ON==1): 
                            
                            response_moment = 1
                            
                            if key_pressed == correct_resp:
                                acc = 1
                                if key_pressed == touchIs[0]: trigg_resp = 101 # correct left button response
                                if key_pressed == touchIs[1]: trigg_resp = 103 # correct right button response 

                            # Incorrect response
                            elif key_pressed != correct_resp:
                                acc = 0
                                if key_pressed == touchIs[0]: trigg_resp = 100 #incorrect left button response
                                if key_pressed == touchIs[1]: trigg_resp = 102 #incorrect right button response 

                            rt = theseKeys.rt 
                            if activate_port is True: port.write(bytes([trigg_resp]))
                            win.flip()
                            break                        
                        
                    # 'Escp' to quit
                    if defaultKeyboard.getKeys(keyList=["escape"], waitRelease=False): core.quit()

                    # Draw circle and dots
                    if (Resp_Signal_ON == 0):
                        if (count_frame <= time_stim_delay):
                            fixation_cross_blue.draw() # draw circle during stim and delay period
                            if (count_frame <= stim_pres_trial): 
                                dots.draw()
                                #core.wait(1)

                    
                    # Trigger end of presentation dots
                    if (Resp_Signal_ON == 0) and (Stim_ON == 1) and (count_frame > stim_pres_trial): 
                        if activate_port is True: win.callOnFlip(port.write, bytes([int(1)]))
                        Stim_ON = 0
                   
                    if (Resp_Signal_ON == 0) and (count_frame > time_stim_delay) : # time_stim_delay 120 + 180 = 300 grame
                        fixation_cross_yellow.draw() # RS
                        trigg_stim_delay = int(trigg_stim) + 100
                        
                        print('RS', count_frame, '- time_stim_delay', time_stim_delay)

                        win.callOnFlip(resp.clock.reset) # reset clock at the end of delay
                        if activate_port is True: win.callOnFlip(port.write, bytes([trigg_stim_delay]))       
                        Resp_Signal_ON = 1

                if not kb_delay:
                    # Draw first dots
                    dots.draw()
                    fixation_cross_blue.draw()

                    resp.clearEvents(eventType='keyboard')
                    kb_delay = 1 

                count_frame += 1
                
                
                # --Save coordinate of dot for each frame to capture motion energy (not used)
                #record frame coordinates
                current_frame =  np.zeros_like(xx, dtype=np.bool) ##empty frame
                y = dots._verticesBase.transpose()
                m, n =  np.round(y * ppd + rayon_cercle * ppd -  dotsize / 2).astype(np.int)
                current_frame[m, n] = 1 ##fill the mesh
                #current_frame = maximum_filter(current_frame,  dotsize) ##build dots of size=dotsize / commented because framedrops
                sparse_current_frame =  sparse.csr_matrix(current_frame)
                dots_coord.append(sparse_current_frame)   


                win.flip()  # function resp.clock.reset synchronized with the frame flip with callOnFlip

        win.recordFrameIntervals = False
        win.flip()  
        dropped[trial_tot-1] = win.nDroppedFrames

        # save dots coordinate for motion energy
        dots_coord = np.array(dots_coord)
        dots_coord = np.stack(dots_coord,axis = -1)

        os.chdir(folder_dots)
        np.save("s"+str(subject_id)+"_trial"+ str(trial_tot) + ".npy",dots_coord)


        # Save data
        data_array = [  subject_id,
                        trial_tot,
                        block_counter,
                        type_block,
                        trial_counter, 
                        coherence_trial,
                        stim_pres_trial,
                        direction_trial,
                        correct_resp,
                        key_pressed,
                        response_moment,
                        acc, # 1 if correct   0 if incorrect
                        rt,
                        trigg_stim,
                        trigg_resp,
                        dropped[trial_tot-1]
                        ]   
                              
        data_container.append(data_array)
        with open(running_filename, 'a') as f_handle: f_handle.write( ','.join([str(data) for data in data_array])+'\n' ) # "with open" function close automatically

        trial_counter += 1
        trial_tot += 1

    # Pause instruction at the end of block
    if block_counter < len(general_seq): 
        # save between block
        data_container_tmp = np.array(data_container)
        np.save(folder+filename, data_container_tmp) 
        # draw instructions
        instruction_exp(instructions_pause, touchIs)
   
    block_counter += 1
    # End du bloc --------------------------------

os.chdir(folder)
data_container = np.array(data_container)
np.save(folder+filename, data_container)    



np.save("s"+str(subject_id)+ "_dropped_frames.npy",dropped)

# End
win.flip()
fin.draw()
win.flip()
core.wait(4)

# Shutting down
win.close()
core.quit()


# [X]