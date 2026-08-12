#!/usr/bin/env python3
"""Signature sweep v2: for harvested (poly -> 24Tt, r_known), find groups with
OPEN r values; sweep f(x)+c across critical-value intervals so r hits the open
values (Hilbert irreducibility: group usually preserved for generic c when the
base poly is generic for its group).

Input:
  harvest.json  [{"poly": "a0,...,a24", "t": N, "r": N, "label": "24T..."}]
  targets.tsv   t<TAB>r<TAB>label
Output: sweep candidates to stdout.
"""
import json, sys, random
import numpy as np
from cypari2 import Pari

pari = Pari()
X = pari("x")

def load_harvest(path):
    out = []
    for h in json.load(open(path)):
        c = [int(z) for z in h["poly"].split(",")]
        out.append((c, h["t"], h["r"]))
    return out

def load_targets(path):
    by_t = {}
    with open(path) as f:
        next(f)
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) < 2: continue
            by_t.setdefault(int(parts[0]), set()).add(int(parts[1]))
    return by_t

def critical_values(coeffs):
    """Values of f at real critical points (roots of f')."""
    d = len(coeffs) - 1
    der = np.polyder(coeffs[::-1])  # numpy descending order
    crit = np.roots(der)
    crit = crit[np.abs(crit.imag) < 1e-6].real
    vals = np.polyval(coeffs[::-1], crit)
    return np.sort(vals)

def sweep_cs(coeffs, n_tail=8):
    """c values to try: midpoints of critical-value intervals + tails."""
    vals = critical_values(coeffs)
    cs = []
    if len(vals) == 0:
        return [0]
    cs.append(vals[0] - (vals[-1] - vals[0] + 10))
    for i in range(len(vals) - 1):
        cs.append((vals[i] + vals[i + 1]) / 2)
    cs.append(vals[-1] + (vals[-1] - vals[0] + 10))
    # spread tails further
    span = vals[-1] - vals[0] + 10
    for i in range(1, n_tail):
        cs.append(vals[0] - span * 2 ** i)
        cs.append(vals[-1] + span * 2 ** i)
    return cs

def count_real_roots(g):
    try:
        return int(pari.polsturm(g))
    except Exception:
        return -1

def is_irreducible(g):
    try:
        return bool(g.polisirreducible())
    except Exception:
        return False

def main():
    rng = random.Random(20260815)
    harvest = load_harvest("data/harvest.json")
    targets = load_targets("data/targets.tsv")
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    by_t = {}
    for c, t, r in harvest:
        by_t.setdefault(t, []).append((c, r))
    gaps = [t for t, rs in targets.items() if t in by_t]
    print(f"harvest groups: {len(by_t)}, with open-r gaps: {len(gaps)}", file=sys.stderr)
    out = []
    seen = set()
    for t in sorted(gaps):
        open_rs = sorted(targets[t])
        for c0, rk in by_t[t]:
            cs = sweep_cs(c0)
            for c in cs:
                cc = c0.copy()
                cc[0] = c0[0] + int(round(c))
                if cc[0] == 0: continue
                key = tuple(cc)
                if key in seen: continue
                g = pari.Polrev(cc)
                if not is_irreducible(g): continue
                r = count_real_roots(g)
                if r in open_rs:
                    seen.add(key)
                    out.append((cc, t, r))
                    break  # one candidate per target r hit per base poly
        print(f"  t={t} open_r={open_rs}: {sum(1 for o in out if o[1]==t)} cands", file=sys.stderr)
    rng.shuffle(out)
    out = out[:n_max]
    print(f"total: {len(out)}", file=sys.stderr)
    for c, t, r in out:
        print(",".join(str(z) for z in c) + f" # sweep t={t} r={r}")

if __name__ == "__main__":
    main()
