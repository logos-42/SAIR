#!/usr/bin/env python3
"""Merge + dedupe all unsubmitted lift candidates into batch15 (≤ 1000)."""
import sys

files = [
    "data/liftE_c13_full.txt", "data/liftE_c13_full_v2.txt", "data/liftE_c13_full_v2b.txt",
    "data/liftE_c13_mix.txt", "data/liftE_c13_neg.txt", "data/liftE_c9xc2_full.txt",
    "data/liftE_s3x2_full.txt", "data/liftE_s3x2_full_v2b.txt",
    "data/liftB_b1_full.txt", "data/liftB_b1_mix.txt", "data/liftB_b1_neg.txt",
    "data/liftB_b2_full_n39.txt", "data/liftB_b2_full_n45.txt", "data/liftB_b2_full_n56.txt",
]
seen = set()
out = []
for f in files:
    try:
        lines = [l.strip() for l in open(f) if l.strip()]
    except FileNotFoundError:
        continue
    for l in lines:
        poly = l.split("#")[0].strip()
        if poly in seen:
            continue
        seen.add(poly)
        out.append(l)
with open("data/explore_batch15.txt", "w") as fh:
    fh.write("\n".join(out) + "\n")
print(f"merged {len(out)} unique candidates")
