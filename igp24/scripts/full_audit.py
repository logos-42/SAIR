#!/usr/bin/env python3
"""Full audit: every (t, r) ever submitted (all batch logs) vs open targets.

Scans /tmp/submit_*.log and data/*_submit.log for idx lines with labels;
matches against new_targets.json.  This is the definitive check for any
missed open hit across ALL families (helix/lift/sphere/sparse/compositum).
"""
import json, re, glob, os, collections

nt = json.load(open("data/new_targets.json"))
by_t = collections.defaultdict(set)
for x in nt:
    by_t[x["t"]].add(x["r"])

logs = sorted(glob.glob("/tmp/submit_*.log")) + sorted(glob.glob("data/*submit*.log"))
all_pairs = collections.Counter()
hits = collections.Counter()
seen_logs = 0
for path in logs:
    try:
        txt = open(path, errors="ignore").read()
    except Exception:
        continue
    if "idx=" not in txt and "24T" not in txt:
        continue
    seen_logs += 1
    for m in re.finditer(r"(?:idx=\d+\s+)?(24T\d+)\s+r=(\d+)", txt):
        t, r = m.group(1), int(m.group(2))
        all_pairs[(t, r)] += 1
        if t in by_t and r in by_t[t]:
            hits[(t, r)] += 1

print(f"scanned {seen_logs} logs, {len(all_pairs)} distinct (t,r), "
      f"{sum(all_pairs.values())} total submissions")
print("OPEN HITS:", dict(hits) if hits else "NONE")
# distribution of all submitted labels vs open labels
submitted_t = set(t for t, r in all_pairs)
open_t = set(by_t.keys())
overlap = submitted_t & open_t
print(f"submitted {len(submitted_t)} distinct groups; {len(open_t)} open groups; "
      f"overlap {len(overlap)}: {sorted(overlap)[:20]}")
