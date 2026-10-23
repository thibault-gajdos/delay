"""shardio.py -- read sweep records from shard files.

Every record of a draw depends only on its index and on the leaks, never on how the run
was split or interrupted, so two records of the same draw in the same directory are
identical: they can only arise from two launches running at once.  load() keeps one of
each and skips an incomplete last line, so a directory left by a crash or by a duplicate
launch is still read correctly.
"""
import glob, json, os


def load(pattern):
    """All records in the files matching `pattern`, one per draw, sorted by draw."""
    out = {}
    for path in sorted(glob.glob(pattern)):
        for line in open(path):
            if not line.endswith('\n'):
                continue                       # incomplete last line of a crashed run
            r = json.loads(line)
            out.setdefault(r['draw'], r)
    return [out[d] for d in sorted(out)]
