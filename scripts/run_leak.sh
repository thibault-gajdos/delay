#!/bin/bash
# Leak-robustness reruns (reviewer C2.3), run natively: the whole grid in one go.
#
#   ./run_leak.sh
#
# Same grid, draws, trial count, shard count and output layout as chunk_leak.sh, so the
# two can be used on the same output directory: gamma_x, gamma_q in {0, 0.25, 0.5, 1} /s,
# every combination but (0, 0), each re-running the draws admissible in either model at
# zero leak (results/leak/adm_idx.json, from leak_index.py).
set -eu
cd "$(dirname "$0")"
NDRAW=1200; NTR=10000; P=4; G="0 0.25 0.5 1"
IDX=results/leak/adm_idx.json
[ -f "$IDX" ] || python3 leak_index.py
for gx in $G; do for gq in $G; do
  [ "$gx" = 0 ] && [ "$gq" = 0 ] && continue
  D=results/leak/gx${gx}_gq${gq}
  mkdir -p "$D"
  for i in $(seq 0 $((P-1))); do
    IDX=$IDX SHARD_N=$P SHARD_I=$i python3 fused_sweep.py $NDRAW $NTR \
      "$D/shard_$i.jsonl" $gx $gq >> "$D/log_$i.txt" 2>&1 &
  done
  wait
  echo "gamma_x $gx, gamma_q $gq: $(cat "$D"/shard_*.jsonl | wc -l | tr -d ' ') draws"
done; done
