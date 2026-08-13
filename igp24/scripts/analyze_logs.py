#!/usr/bin/env python3
"""Analyze verified batch logs: distinct (t, r) counts."""
import re, sys, collections

for path in sys.argv[1:]:
    counts = collections.Counter()
    try:
        txt = open(path).read()
    except Exception as e:
        print(path, "ERR", e)
        continue
    for m in re.finditer(r"idx=\d+ (24T\d+) r=(\d+)", txt):
        counts[(m.group(1), int(m.group(2)))] += 1
    if not counts:
        print(path, "no hits found")
        continue
    top = counts.most_common(14)
    print(f"=== {path} ({sum(counts.values())} polys, {len(counts)} distinct) ===")
    for (t, r), n in top:
        print(f"  {t} r={r}: {n}")
