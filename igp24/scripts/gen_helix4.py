#!/usr/bin/env python3
"""Helix v4: shape 2·2·3·2 — ∛ layer after TWO rotation layers, √ layer on top.

K0 = Q(ζ3)          (rotation, degree 2)
K1 = K0(√d)         (rotation, degree 4)
K2 = K1(∛a)         (Kummer, degree 12; 3-part via K1's conjugates)
K3 = K2(√b)         (quadratic on top; 2-part amplification via degree-12 base)

Different layer order than v3 → different 2-part/3-part tradeoff.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260825

def is_irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False

def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1

def prim_part(c):
    gg = 0
    for z in c:
        gg = math.gcd(gg, abs(z))
    if gg > 1:
        return [z // gg for z in c]
    return c

def layer(f, bcoeffs, kth, do_reduce=True):
    d = f.poldegree()
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs, Y)
    coeffs = [Z**kth - bcoeffs[0]] + [-bcoeffs[i] for i in range(1, d)]
    g = pari.Pol(coeffs, Y)
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    d2 = kth * d
    if R.poldegree() != d2:
        return None
    c = prim_part([int(z) for z in pari.Vecrev(pari.Pol(R, X))])
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    if do_reduce and d2 <= 12:
        try:
            Pr = pari.polredabs(P)
        except Exception:
            Pr = P
        if Pr.poldegree() == d2:
            P = Pr
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return cc

def rand_elem(f, rng, bound=4):
    d = f.poldegree()
    for _ in range(40):
        bc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    while made < n_target and tries < 3000:
        tries += 1
        f0 = pari.Polrev([1, 1, 1])  # Q(ζ3)
        bc = rand_elem(f0, rng)
        if bc is None: continue
        c = layer(f0, bc, 2)
        if c is None: continue
        f1 = pari.Polrev(c)  # degree 4
        bc = rand_elem(f1, rng, bound=3)
        if bc is None: continue
        c = layer(f1, bc, 3)  # ∛ on degree-4 base
        if c is None: continue
        f2 = pari.Polrev(c)  # degree 12
        bc = rand_elem(f2, rng, bound=3)
        if bc is None: continue
        c = layer(f2, bc, 2, do_reduce=False)  # √ on top
        if c is None: continue
        if len(c) != 25: continue
        lc = c[24]
        if lc == 1: pass
        elif lc == -1: c = [-z for z in c]
        elif all(z % lc == 0 for z in c): c = [z // lc for z in c]
        else:
            o2 = []
            for k in range(24):
                o2.append(c[k] * lc ** (23 - k))
            o2.append(1)
            c = o2
        key = tuple(c)
        if key in seen: continue
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        r = real_count(g)
        if r < 0 or r % 2 != 0: continue
        seen.add(key)
        out.append((c, "helix4"))
        made += 1
    print(f"helix4: made={made} tries={tries}", file=sys.stderr)
    for c, tag in out:
        r = real_count(pari.Polrev(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
