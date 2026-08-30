#!/bin/bash
# Waits for the behavioral G^2 calibration to finish (both models, run separately via
# supervise_calib.sh), then refreshes the full EMG+behavioral calibration
# (results/fit_model1.json / results/fit_model2.json) under the 1024 Hz DT.
# calibrate.py checkpoints its own progress in results/fit_model{N}.json, so this is
# safe to interrupt and resume.
cd "$HOME/thib/projects/delay/scripts"
export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"

while [ ! -f results/calibq_m1.json ] || [ ! -f results/calibq_m2.json ]; do
  sleep 20
done
echo "$(date): behavioral G2 calibration done, starting EMG recalibration under 1024Hz" >> results/logs/emg_recal.log

for chunk in 1 2 3 4; do
  [ -f results/DONE_EMG_RECAL ] && break
  uv run python3 calibrate.py 90 1 >> results/logs/emg_recal_m1.log 2>&1
  uv run python3 calibrate.py 90 2 >> results/logs/emg_recal_m2.log 2>&1
  echo "$(date): chunk $chunk done" >> results/logs/emg_recal.log
done
touch results/DONE_EMG_RECAL
echo "$(date): EMG recalibration finished" >> results/logs/emg_recal.log
