#!/bin/bash
# One bounded slice of the parameter-space sweep, for a shell that gets a fresh process
# per call with a hard time limit (no long background job possible).
#
#   ./chunk_sweep.sh 165        # run for at most 165 s on 4 shards, then stop
#
# fused_sweep.py flushes one line per draw and resumes by line count, so calling this
# repeatedly walks the sweep to completion; a draw interrupted by the time limit is
# simply recomputed on the next call.  The shard count is fixed at 4 so that the
# shard -> draw mapping does not depend on the machine.  On a machine that can run the
# job in one go, use run_sweep.sh instead (it uses its own shard count; do not mix the
# two on the same output directory).
set -u
cd "$(dirname "$0")"
SECS=${1:-165}
NDRAW=2400; NTR=10000; P=4; OUTDIR=results/shards
mkdir -p "$OUTDIR"
echo "$P" > "$OUTDIR/.shards"
for i in $(seq 0 $((P-1))); do
  SHARD_N=$P SHARD_I=$i timeout "$SECS" python3 fused_sweep.py $NDRAW $NTR \
    "$OUTDIR/shard_$i.jsonl" 0.0 >> "$OUTDIR/log_$i.txt" 2>&1 &
done
wait
echo "$(cat "$OUTDIR"/shard_*.jsonl 2>/dev/null | wc -l | tr -d ' ') / $NDRAW draws"
