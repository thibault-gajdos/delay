#!/bin/bash
# One bounded slice of the leak rerun, for environments where a long-running background
# job is not possible (the sandboxed shell this repository is sometimes driven from gets a
# fresh process per call, with a hard time limit).
#
#   ./chunk_leak.sh 0.25 150       # run gamma = 0.25 for at most 150 s across all shards
#
# leak_sweep.py flushes one line per draw and resumes by line count, so calling this
# repeatedly with the same arguments walks the draw list to completion. Prints progress.
# On a machine that can run the job to completion in one go, use run_leak.sh instead.
set -u
cd "$(dirname "$0")"
G=${1:-0.25}; SECS=${2:-150}
tag=${G/./}
OUTDIR=results/leak
P=$(nproc 2>/dev/null || sysctl -n hw.physicalcpu 2>/dev/null || echo 4)
mkdir -p "$OUTDIR"
[ -f "$OUTDIR/adm_idx.json" ] || python3 leak_index.py
N=$(python3 -c "import json;print(len(json.load(open('$OUTDIR/adm_idx.json'))))")
for i in $(seq 0 $((P-1))); do
  GAMMA=$G OUT="$OUTDIR/leak_g${tag}_sh$i.jsonl" SHARD_N=$P SHARD_I=$i \
    timeout "$SECS" python3 leak_sweep.py >> "$OUTDIR/log_g${tag}_$i.txt" 2>&1 &
done
wait
DONE=$(cat "$OUTDIR"/leak_g${tag}_sh*.jsonl 2>/dev/null | wc -l | tr -d ' ')
echo "gamma $G: $DONE / $N draws"
if [ "$DONE" -ge "$N" ]; then
  cat "$OUTDIR"/leak_g${tag}_sh*.jsonl > "$OUTDIR/leak_g${tag}.jsonl"
  echo "gamma $G COMPLETE -> $OUTDIR/leak_g${tag}.jsonl"
fi
