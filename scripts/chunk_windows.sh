#!/bin/bash
# One bounded slice of augment_windows.py, for a shell with a hard time limit per call.
#   ./chunk_windows.sh 100        # run for at most 100 s across all shards
# Resumable: call repeatedly until it reports the full draw count.
set -u
cd "$(dirname "$0")"
SECS=${1:-100}
OUTDIR=results/windows
P=$(nproc 2>/dev/null || sysctl -n hw.physicalcpu 2>/dev/null || echo 4)
mkdir -p "$OUTDIR"
N=$(python3 -c "import json;print(len(json.load(open('results/leak/adm_idx.json'))))")
for i in $(seq 0 $((P-1))); do
  SHARD_N=$P SHARD_I=$i timeout "$SECS" python3 augment_windows.py "$OUTDIR/sh$i.jsonl" \
    >> "$OUTDIR/log_$i.txt" 2>&1 &
done
wait
DONE=$(cat "$OUTDIR"/sh*.jsonl 2>/dev/null | wc -l | tr -d ' ')
echo "$DONE / $N draws"
if [ "$DONE" -ge "$N" ]; then
  cat "$OUTDIR"/sh*.jsonl > results/window_effects.jsonl
  echo "COMPLETE -> results/window_effects.jsonl"
fi
