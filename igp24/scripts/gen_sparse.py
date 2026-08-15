#!/usr/bin/env python3
"""Sparse polynomial scanner: x^24 + Σ c_i x^i (few terms, small coeffs).

Rationale: tower products are dense (huge closure groups, all covered).
Sparse polys (2-5 non-zero terms, coeffs ±1,±2) have constrained Galois
groups via Newton-polygon / ramification arguments — a family never
systematically scanned before.  r controlled by Sturm.

Families:
  F1  x^24 + a·x^k + b            (trinomial, k ∈ 1..23)
  F2  x^24 + a·x^k + b·x^m + c    (4-term)
  F3  x^24 + a·x^12 + b           (even/odd symmetry → r structure)
  F4  x^24 + a·x^k + b·x^(24−k) + c  (palindromic → reciprocal roots)
"""
import random, sys, math
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260901
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 60
FAMILY = int(sys.argv[3]) if len(sys.argv) > 3 else 0  # 0=all, 1=F1, 2=F2, 3=F3, 4=F4

def is_irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False

def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    fams = [FAMILY] if FAMILY else (1, 2, 3, 4)
    while made < N_TARGET and tries < 300 * N_TARGET:
        tries += 1
        fam = rng.choice(fams)
        c = [0] * 25
        c[24] = 1
        if fam == 1:
            k = rng.randint(1, 23)
            a = rng.choice((1, -1, 2, -2, 3, -3))
            b = rng.choice((1, -1, 2, -2, 3, -3))
            c[k] = a
            c[0] = b
        elif fam == 2:
            ks = rng.sample(range(1, 24), 2)
            for k in ks:
                c[k] = rng.choice((1, -1, 2, -2))
            c[0] = rng.choice((1, -1, 2, -2))
        elif fam == 3:
            a = rng.choice((1, -1, 2, -2, 3, -3))
            b = rng.choice((1, -1, 2, -2, 3, -3))
            c[12] = a
            c[0] = b
        else:  # F4 palindromic
            k = rng.randint(1, 11)
            a = rng.choice((1, -1, 2, -2, 3, -3))
            b = rng.choice((1, -1, 2, -2, 3, -3))
            c[k] = a
            c[24 - k] = a
            c[0] = b
        if c[0] == 0:
            continue
        key = tuple(c)
        if key in seen:
            continue
        g = pari.Polrev(c)
        if not is_irred(g):
            continue
        r = real_count(g)
        if r < 0 or r % 2 != 0:
            continue
        seen.add(key)
        out.append((c, r, fam))
        made += 1
    print(f"sparse: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    fdist = {}
    for c, r, fam in out:
        rdist[r] = rdist.get(r, 0) + 1
        fdist[fam] = fdist.get(fam, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    print(f"family分布: {fdist}", file=sys.stderr)
    for c, r, fam in out:
        print(",".join(str(z) for z in c) + f" # sparseF{fam} r={r}")

if __name__ == "__main__":
    main()
