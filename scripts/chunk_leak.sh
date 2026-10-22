#!/bin/bash
# One bounded slice of the leak reruns (reviewer C2.3), for a shell that gets a fresh
# process per call with a hard time limit.
#
#   ./chunk_leak.sh 165        # work for at most 165 s on 4 shards, then stop
#
# The leaks of the two accumulators are varied independently over a full grid,
# gamma_x, gamma_q in {0, 0.25, 0.5, 1} /s, every combination but (0, 0), which is the
# main sweep. Each combination re-runs the draws admissible in either model at zero leak
# (results/leak/adm_idx.json, from leak_index.py), at the main sweep's trial count and
# seeds, into results/leak/gx<gamma_x>_gq<gamma_q>/. Each call works on the first
# combination not yet complete; fused_sweep.py resumes by line count, so calling this
# repeatedly walks the whole grid. On a machine that can run the job in one go, use
# run_leak.sh instead.
set -u
cd "$(dirname "$0")"
SECS=${1:-165}
NDRAW=1200; NTR=10000; P=4; G="0 0.25 0.5 1"
IDX=results/leak/adm_idx.json
[ -f "$IDX" ] || python3 leak_index.py
N=$(python3 -c "import json;print(len(json.load(open('$IDX'))))")
for gx in $G; do for gq in $G; do
  [ "$gx" = 0 ] && [ "$gq" = 0 ] && continue
  D=results/leak/gx${gx}_gq${gq}
  DONE=$(cat "$D"/shard_*.jsonl 2>/dev/null | wc -l | tr -d ' ')
  [ "$DONE" -ge "$N" ] && continue
  mkdir -p "$D"
  for i in $(seq 0 $((P-1))); do
    IDX=$IDX SHARD_N=$P SHARD_I=$i timeout "$SECS" python3 fused_sweep.py $NDRAW $NTR \
      "$D/shard_$i.jsonl" $gx $gq >> "$D/log_$i.txt" 2>&1 &
  done
  wait
  echo "gamma_x $gx, gamma_q $gq: $(cat "$D"/shard_*.jsonl | wc -l | tr -d ' ') / $N draws"
  exit 0
done; done
echo "leak grid complete"
