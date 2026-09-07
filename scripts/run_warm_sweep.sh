#!/bin/bash
# Stationary-start parameter sweep, run natively (macOS or Linux).
#
#   ./run_warm_sweep.sh          # uses all physical cores
#   ./run_warm_sweep.sh 8        # 8 shards
#
# fused_sweep.py appends one flushed JSON line per draw and resumes by counting the
# lines already in its shard file, so the run can be interrupted and restarted freely.
# Changing the shard count invalidates the shard->draw mapping, so the script wipes
# previous output if the count differs from the recorded one.
set -eu
cd "$(dirname "$0")"
NDRAW=1200; NTR=4000; OUTDIR=results/warm
if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then
  DEF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
else
  DEF=$(nproc 2>/dev/null || echo 4)
fi
P=${1:-$DEF}
mkdir -p "$OUTDIR"
PREV=$(cat "$OUTDIR/.shards" 2>/dev/null || echo "$P")
if [ "$PREV" != "$P" ]; then
  echo "shard count changed ($PREV -> $P): wiping $OUTDIR and starting fresh"
  rm -f "$OUTDIR"/shard_*.jsonl "$OUTDIR"/log_*.txt
fi
echo "$P" > "$OUTDIR/.shards"
echo "launching $P shards over $NDRAW draws x $NTR trials/condition"
pids=()
for i in $(seq 0 $((P-1))); do
  SHARD_N=$P SHARD_I=$i python3 fused_sweep.py $NDRAW $NTR "$OUTDIR/shard_$i.jsonl" 0.0 \
      > "$OUTDIR/log_$i.txt" 2>&1 &
  pids+=($!)
done
echo "pids: ${pids[*]}"
echo "progress:  cat $OUTDIR/shard_*.jsonl | wc -l    (target $NDRAW)"
wait
echo "done: $(cat "$OUTDIR"/shard_*.jsonl | wc -l) draws"
