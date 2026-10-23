"""repair_shards.py -- make shard files safe to resume from after a crash.

  python3 repair_shards.py FILE...

A run killed in the middle of a write can leave an incomplete last line.  This drops it,
in every file given, so that the draw is simply recomputed on the next run.  It is called
by run_sweep.sh and run_leak.sh before every launch, on every shard file of the output
directory, including files no process of the new launch will open (when the shard count
has changed).
"""
import os, sys


def repair(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    with open(path, 'rb+') as f:
        data = f.read()
        if data.endswith(b'\n'):
            return False
        f.seek(0)
        f.truncate(data.rfind(b'\n') + 1)
        return True


if __name__ == '__main__':
    for p in sys.argv[1:]:
        if repair(p):
            print(f"repaired {p}: dropped an incomplete last line")
