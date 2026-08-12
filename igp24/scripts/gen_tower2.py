#!/usr/bin/env python3
"""Tower generator v2: nested quadratic tower over an S3 cubic field.

L = K3(√u)(√v)(√w), K3 = Q(α) cubic (S3 or cyclic), u,v,w ∈ intermediate fields.
Each layer: minimal poly via resultant Res_y(f(y), z²−b(y)) — degree doubles.
Nested √-tower gives non-abelian 2-groups (D4/Q8/… ⋊ S3 data), r controlled by
signs of σ(b) per layer (can reach r=24 with all-positive choices).

Shapes:
  n2:  3 layers of √ over cubic     -> degree 3·2·2·2 = 24
  n1:  cubic + √ + quartic-ish?     (not used)
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260820

def is_irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False

def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1

def monicize(c):
    lc = c[24]
    if lc == 1: return c
    if lc == -1: return [-z for z in c]
    if all(z % lc == 0 for z in c): return [z // lc for z in c]
    out = []
    for k in range(24):
        out.append(c[k] * lc ** (23 - k))
    out.append(1)
    return out

def sq_layer(f, bcoeffs):
    """Res_y(f(y), z²−b(y)) with z→x; monic degree-2d poly (if irreducible)."""
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs, Y)
    coeffs = [Z**2 - bcoeffs[0]] + [-bcoeffs[i] for i in range(1, len(bcoeffs))]
    g = pari.Pol(coeffs, Y)
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    d = 2 * (f.poldegree())
    if R.poldegree() != d:
        return None
    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
    gg = 0
    for z in c:
        gg = math.gcd(gg, abs(z))
    if gg > 1:
        c = [z // gg for z in c]
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    # reduce to small-coefficient defining poly of the same field
    try:
        Pr = pari.polredabs(P)
    except Exception:
        Pr = P
    if Pr.poldegree() != d:
        Pr = P
    cc = [int(z) for z in pari.Vecrev(pari.Pol(Pr, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return cc

def random_b(f, rng, bound=6):
    """Random element of Q[x]/(f) as coeff list (degree < deg f)."""
    d = f.poldegree()
    for _ in range(20):
        bc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def cubic(rng, bound=9):
    """Random irreducible cubic (S3 typically; r=1 or 3)."""
    for _ in range(200):
        c = [rng.randint(-bound, bound), rng.randint(-bound, bound), 0, 1]
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if is_irred(g):
            return c
    return None

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    while made < n_target and tries < 3000:
        tries += 1
        c3 = cubic(rng)
        if c3 is None: continue
        f = pari.Polrev(c3)  # degree 3
        # layer 1: √u, u ∈ K3
        bc = random_b(f, rng)
        if bc is None: continue
        c = sq_layer(f, bc)
        if c is None: continue
        f2 = pari.Polrev(c)  # degree 6
        # layer 2: √v, v ∈ K6
        bc = random_b(f2, rng, bound=5)
        if bc is None: continue
        c = sq_layer(f2, bc)
        if c is None: continue
        f3 = pari.Polrev(c)  # degree 12
        # layer 3: √w, w ∈ K12
        bc = random_b(f3, rng, bound=4)
        if bc is None: continue
        c = sq_layer(f3, bc)
        if c is None: continue
        if len(c) != 25: continue
        c = monicize(c)
        key = tuple(c)
        if key in seen: continue
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        r = real_count(g)
        if r < 0 or r % 2 != 0: continue
        seen.add(key)
        out.append((c, f"tower2 cube={tuple(c3)}"))
        made += 1
    print(f"tower2: made={made} tries={tries}", file=sys.stderr)
    for c, tag in out:
        r = real_count(pari.Polrev(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
