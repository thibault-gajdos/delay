#!/bin/bash
# Leak-robustness reruns (reviewer C2.3), run natively, one shard per core.
#
#   ./run_leak.sh          # one shard per performance core
#   ./run_leak.sh 8        # 8 shards
#
# gamma_x, gamma_q in {0, 0.25, 0.5, 1} /s, every combination but (0, 0), which is the main
# sweep. Each combination runs all the draws of the main sweep, into
# results/leak/gx<gamma_x>_gq<gamma_q>/, and admissibility is evaluated at that
# combination. Draws already recorded are skipped, so the run can be interrupted and
# restarted, with any shard count. It does not depend on the main sweep and can run
# before, after or alongside it.
set -eu
cd "$(dirname "$0")"
NDRAW=1200; NTR=40000; G="0 0.25 0.5 1"
if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then
  DEF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
else
  DEF=$(nproc 2>/dev/null || echo 4)
fi
P=${1:-$DEF}
for gx in $G; do for gq in $G; do
  [ "$gx" = 0 ] && [ "$gq" = 0 ] && continue
  D=results/leak/gx${gx}_gq${gq}
  mkdir -p "$D"
  python3 repair_shards.py "$D"/shard_*.jsonl 2>/dev/null || true
  for i in $(seq 0 $((P-1))); do
    DONE="$D/shard_*.jsonl" SHARD_N=$P SHARD_I=$i \
      python3 fused_sweep.py $NDRAW $NTR "$D/shard_$i.jsonl" $gx $gq \
      >> "$D/log_$i.txt" 2>&1 &
  done
  wait
  echo "gamma_x $gx, gamma_q $gq: $(cat "$D"/shard_*.jsonl | wc -l | tr -d ' ') / $NDRAW draws"
done; done
