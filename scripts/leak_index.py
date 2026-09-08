"""leak_index.py -- the draw list the leak reruns walk.

Every draw admissible in EITHER model at gamma = 0, taken from the sweep dump given on
the command line (default: the stationary-start dump).  Written next to the dump's own
leak directory so the cold-start list is not overwritten.

Usage:  python3 leak_index.py [DUMP.jsonl]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = sys.argv[1] if len(sys.argv) > 1 else 'results/sweep_screened_warm.jsonl'
SFX = 'warm' if 'warm' in os.path.basename(DUMP) else 'cold'
OUT = os.path.join(HERE, f'results/leak/{SFX}/adm_idx.json')

recs = [json.loads(l) for l in open(os.path.join(HERE, DUMP))]
idx = sorted(r['draw'] for r in recs
             if r['model1']['admissible'] or r['model2']['admissible'])
n1 = sum(1 for r in recs if r['model1']['admissible'])
n2 = sum(1 for r in recs if r['model2']['admissible'])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(idx, open(OUT, 'w'))
print(f"{DUMP}: {n1} admissible in Model 1, {n2} in Model 2, "
      f"{len(idx)} in either -> {os.path.relpath(OUT, HERE)}")
