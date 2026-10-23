#!/bin/bash
# One bounded slice of the parameter-space sweep, for a shell that gets a fresh process
# per call with a hard time limit (no long background job possible).
#
#   ./chunk_sweep.sh 165        # run for at most 165 s on 4 shards, then stop
#
# fused_sweep.py flushes one line per draw and skips recorded draws, so calling this
# repeatedly walks the sweep to completion; a draw interrupted by the time limit is
# simply recomputed on the next call.  Draws already recorded in any shard file
# are skipped, so this and run_sweep.sh can be mixed freely.  On a machine that can run
# the job in one go, use run_sweep.sh instead.
set -u
cd "$(dirname "$0")"
# `timeout` is GNU coreutils; macOS has it only as gtimeout (brew coreutils)
TO=$(command -v timeout || command -v gtimeout || true)
[ -n "$TO" ] && TO="$TO ${1:-165}" || TO=""
SECS=${1:-165}
NDRAW=1200; NTR=10000; P=4; OUTDIR=results/shards
mkdir -p "$OUTDIR"
for i in $(seq 0 $((P-1))); do
  DONE="$OUTDIR/shard_*.jsonl" SHARD_N=$P SHARD_I=$i $TO python3 fused_sweep.py $NDRAW $NTR \
    "$OUTDIR/shard_$i.jsonl" 0.0 >> "$OUTDIR/log_$i.txt" 2>&1 &
done
wait
echo "$(cat "$OUTDIR"/shard_*.jsonl 2>/dev/null | wc -l | tr -d ' ') / $NDRAW draws"
