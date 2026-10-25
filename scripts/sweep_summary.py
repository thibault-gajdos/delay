"""sweep_summary.py -- concatenate the shards of the main sweep.

Run after run_sweep.sh.  Reads results/shards/shard_*.jsonl (one record per draw, whatever
the sharding) and writes them, sorted by draw, to results/sweep_screened.jsonl, the file
every later step reads.

Usage:  python3 sweep_summary.py
"""
import json
from shardio import load

rows = load("results/shards/shard_*.jsonl")
with open("results/sweep_screened.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
n1 = sum(r['model1']['admissible'] for r in rows)
n2 = sum(r['model2']['admissible'] for r in rows)
print(f"{len(rows)} draws -> results/sweep_screened.jsonl; admissible: Model 1 {n1}, "
      f"Model 2 {n2}")
