#!/bin/bash
# Leak-robustness sweep (reviewer C2.3), run natively, one shard per core.
#
#   ./run_leak.sh          # one shard per performance core
#   ./run_leak.sh 8        # 8 shards
#
# 1200 draws of a separate Latin-hypercube cube in which the leaks gamma_x and gamma_q are
# two more sampled dimensions, uniform over [0, 1] /s (fused_sweep.py, argument `sample`),
# into results/leak/. Draws already recorded are skipped, so the run can be interrupted
# and restarted, with any shard count. It does not depend on the main sweep and can run
# before, after or alongside it.
set -eu
cd "$(dirname "$0")"
NDRAW=1200; NTR=40000; OUTDIR=results/leak
if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then
  DEF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
else
  DEF=$(nproc 2>/dev/null || echo 4)
fi
P=${1:-$DEF}
mkdir -p "$OUTDIR"
python3 repair_shards.py "$OUTDIR"/shard_*.jsonl 2>/dev/null || true
echo "launching $P shards over $NDRAW draws x $NTR trials/condition, leaks sampled"
for i in $(seq 0 $((P-1))); do
  DONE="$OUTDIR/shard_*.jsonl" SHARD_N=$P SHARD_I=$i \
    python3 fused_sweep.py $NDRAW $NTR "$OUTDIR/shard_$i.jsonl" sample \
    >> "$OUTDIR/log_$i.txt" 2>&1 &
done
wait
N=$(cat "$OUTDIR"/shard_*.jsonl | wc -l | tr -d ' ')
echo "done: $N / $NDRAW draws"
[ "$N" -ge "$NDRAW" ]
