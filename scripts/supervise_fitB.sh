#!/bin/bash
# Resumable driver for the all-data fit (calib_fitB.py).  Same pattern as
# supervise_calib.sh: relaunch until results/fitB_${tag}.json exists; per-stage progress
# is checkpointed in results/ckptb_${tag}.json so an interrupted run resumes.
# tag is 'm1' or 'm2'.  ~25 min per model.
tag=$1
cd "$HOME/thib/projects/delay/scripts"
export DELAY_DATA="$HOME/thib/projects/delay/data/86cbw-osfstorage-archive"
while [ ! -f results/fitB_${tag}.json ]; do
  uv run python3 calib_fitB.py ${tag} >> results/logs/logb_${tag}_stdout.txt 2>&1
  sleep 2
done
echo "DONE $tag" >> results/logs/logb_${tag}_stdout.txt
