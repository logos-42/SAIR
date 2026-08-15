#!/usr/bin/env python3
"""Helix v7: S3 cubic + MONOMIAL sqrt layers (conjugation-compressed closure).

24T1958 = (D8:(C4×C4)):S3, order 768 = 24·32 — SMALL closure. Random sqrt
elements give independent conjugates → huge closure (2^24). Monomial elements
b = c·α^k (α = defining root) have correlated conjugates → small 2-group.

Shape: K3 (S3 cubic, totally real) + 3 monomial sqrt layers (3·2·2·2 = 24).
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260816

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

def monic_ize(c):
    lc = c[-1]
    if lc == 1: return c
    if lc == -1: return [-z for z in c]
    if all(z % lc == 0 for z in c): return [z // lc for z in c]
    d = len(c) - 1
    out = [c[k] * lc ** (d - 1 - k) for k in range(d)]
    out.append(1)
    return out

def sqrt_layer(f, bcoeffs, do_polred=False):
    d = int(f.poldegree())
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs[::-1], Y)
    g = Z**2 - b
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    if R.poldegree() != 2 * d:
        return None
    c = prim_part([int(z) for z in pari.Vecrev(pari.Pol(R, X))])
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    if do_polred and 2 * d <= 12:
        try:
            Pr = pari.polredabs(P)
            if Pr.poldegree() == 2 * d:
                P = Pr
        except Exception:
            pass
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return monic_ize(cc)

def mono_elem(f, rng, k_max=None, shift=(1, 3)):
    """b = c·α^k + s: single monomial + small shift. Correlated conjugates."""
    d = int(f.poldegree())
    if k_max is None:
        k_max = d
    for _ in range(40):
        k = rng.randint(1, k_max - 1)
        c = rng.choice((1, 2, 3, 5))
        s = rng.randint(shift[0], shift[1])
        bc = [0] * d
        bc[k] = c
        bc[0] += s
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def s3_cubic(rng, bound=7):
    """Random totally-real S3 cubic (disc > 0, non-square)."""
    for _ in range(400):
        c = [rng.randint(-bound, bound), rng.randint(-bound, bound), 0, 1]
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        disc = int(pari.poldisc(g))
        if disc > 0 and real_count(g) == 3:
            # check non-square disc (S3 not C3)
            r = pari.issquare(disc)
            if not r:
                return c
    return None

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    while made < n_target and tries < 500 * n_target:
        tries += 1
        c3 = s3_cubic(rng)
        if c3 is None: continue
        f0 = pari.Polrev(c3)
        bc = mono_elem(f0, rng, k_max=3)
        if bc is None: continue
        c = sqrt_layer(f0, bc, do_polred=True)
        if c is None: continue
        f1 = pari.Polrev(c)
        bc = mono_elem(f1, rng, k_max=6)
        if bc is None: continue
        c = sqrt_layer(f1, bc, do_polred=True)
        if c is None: continue
        f2 = pari.Polrev(c)
        bc = mono_elem(f2, rng, k_max=12)
        if bc is None: continue
        c = sqrt_layer(f2, bc)
        if c is None: continue
        if len(c) != 25: continue
        if c[0] == 0: continue
        key = tuple(c)
        if key in seen: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        r = int(pari.polsturm(g))
        seen.add(key)
        out.append((c, r))
        made += 1
    print(f"helix7: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # helix7 r={r}")

if __name__ == "__main__":
    main()
