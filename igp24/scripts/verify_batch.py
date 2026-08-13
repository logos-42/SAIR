#!/usr/bin/env python3
"""Verify helix batch10 structure: 25 coeffs, monic, a0!=0, irred, r label == polsturm."""
import sys
from cypari2 import Pari
pari = Pari()
path = sys.argv[1] if len(sys.argv) > 1 else "data/explore_batch10.txt"
lines = [l.strip() for l in open(path) if l.strip()]
bad = 0
r_vals = {}
monic = 0
for i, l in enumerate(lines):
    body = l.split("#")[0].strip().split(",")
    if len(body) != 25:
        bad += 1
        continue
    c = [int(z) for z in body]
    if c[24] != 1:
        bad += 1
        continue
    monic += 1
    if c[0] == 0:
        bad += 1
        continue
    g = pari.Polrev(c)
    if not g.polisirreducible():
        bad += 1
        continue
    r_label = int(l.split("r=")[1].split()[0])
    r_act = int(pari.polsturm(g))
    if r_label != r_act:
        bad += 1
        continue
    r_vals[r_label] = r_vals.get(r_label, 0) + 1
print(f"lines={len(lines)} monic={monic} bad={bad} r_dist={dict(sorted(r_vals.items()))}")
print("VERDICT:", "PASS" if bad == 0 and len(lines) > 0 else "FAIL")
