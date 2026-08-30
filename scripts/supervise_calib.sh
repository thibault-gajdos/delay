#!/bin/bash
# Resumable driver for the behavioural G^2 calibration (calib_behav_q.py). Re-launches
# the Python process until results/calibq_${tag}.json exists; calib_behav_q.py itself
# checkpoints per-stage progress in results/ckptq_${tag}.json, so an interrupted run
# (closed terminal, killed process) resumes from the last completed stage rather than
# restarting. tag is 'm1' or 'm2'.
tag=$1
cd "$HOME/mnt/delay/scripts"
while [ ! -f results/calibq_${tag}.json ]; do
  python3 calib_behav_q.py ${tag} >> results/logs/logq_${tag}_stdout.txt 2>&1
  sleep 2
done
echo "DONE $tag" >> results/logs/logq_${tag}_stdout.txt
