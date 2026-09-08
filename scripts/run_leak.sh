#!/bin/bash
# Leak-robustness rerun under the STATIONARY START (reviewer C2.3), run natively.
#
#   ./run_leak.sh          # uses all physical cores
#   ./run_leak.sh 8        # 8 shards
#
# Re-runs, at each of three leak rates, exactly the draws that were admissible in either
# model at zero leak under the stationary start (results/leak/adm_idx.json, written
# by leak_index.py -- 220 draws, against 222 under the cold start). Output goes to
# results/leak/, so the cold-start dumps the response letter currently quotes are
# left alone.
#
# leak_sweep.py appends one flushed JSON line per draw and resumes by counting the lines
# already in its shard file, so the run can be interrupted and restarted freely. Changing
# the shard count invalidates the shard->draw mapping, so the script wipes previous output
# if the count differs from the recorded one.
#
# Cost: three leak rates x 220 draws x 2 models x 6 conditions x 4000 trials, i.e. about
# 55% of a full 1200-draw sweep per leak rate.
set -eu
cd "$(dirname "$0")"
OUTDIR=results/leak
if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then
  DEF=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
else
  DEF=$(nproc 2>/dev/null || echo 4)
fi
P=${1:-$DEF}
mkdir -p "$OUTDIR"
[ -f "$OUTDIR/adm_idx.json" ] || python3 leak_index.py
PREV=$(cat "$OUTDIR/.shards" 2>/dev/null || echo "$P")
if [ "$PREV" != "$P" ]; then
  echo "shard count changed ($PREV -> $P): wiping shard files in $OUTDIR"
  rm -f "$OUTDIR"/leak_g*_sh*.jsonl "$OUTDIR"/log_*.txt
fi
echo "$P" > "$OUTDIR/.shards"
NDRAW=$(python3 -c "import json;print(len(json.load(open('$OUTDIR/adm_idx.json'))))")
echo "launching $P shards over $NDRAW draws, at gamma = 0.25, 0.50, 1.00"
for G in 0.25 0.50 1.00; do
  tag=${G/./}
  echo "--- gamma $G ---"
  pids=()
  for i in $(seq 0 $((P-1))); do
    GAMMA=$G OUT="$OUTDIR/leak_g${tag}_sh$i.jsonl" SHARD_N=$P SHARD_I=$i \
      python3 leak_sweep.py > "$OUTDIR/log_g${tag}_$i.txt" 2>&1 &
    pids+=($!)
  done
  echo "pids: ${pids[*]}   progress:  cat $OUTDIR/leak_g${tag}_sh*.jsonl | wc -l   (target $NDRAW)"
  wait
  cat "$OUTDIR"/leak_g${tag}_sh*.jsonl > "$OUTDIR/leak_g${tag}.jsonl"
  echo "gamma $G done: $(wc -l < "$OUTDIR/leak_g${tag}.jsonl") draws"
done
echo
echo "all done. now:  python3 leak_summary.py"
