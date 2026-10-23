#!/bin/bash
# The parameter-space sweep, run natively (macOS or Linux), one shard per core.
#
#   ./run_sweep.sh          # one shard per performance core
#   ./run_sweep.sh 8        # 8 shards
#
# fused_sweep.py appends one flushed JSON line per draw and skips every draw already
# recorded in any shard file of OUTDIR, so the run can be interrupted, restarted, and
# restarted with a different shard count, without losing or duplicating a draw.
set -eu
cd "$(dirname "$0")"
NDRAW=1200; NTR=10000; OUTDIR=results/shards
if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then
  DEF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
else
  DEF=$(nproc 2>/dev/null || echo 4)
fi
P=${1:-$DEF}
mkdir -p "$OUTDIR"
echo "launching $P shards over $NDRAW draws x $NTR trials/condition"
for i in $(seq 0 $((P-1))); do
  DONE="$OUTDIR/shard_*.jsonl" SHARD_N=$P SHARD_I=$i \
    python3 fused_sweep.py $NDRAW $NTR "$OUTDIR/shard_$i.jsonl" 0 0 \
    >> "$OUTDIR/log_$i.txt" 2>&1 &
done
wait
N=$(cat "$OUTDIR"/shard_*.jsonl | wc -l | tr -d ' ')
echo "done: $N / $NDRAW draws"
[ "$N" -ge "$NDRAW" ]
