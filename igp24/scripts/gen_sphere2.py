#!/usr/bin/env python3
"""Sphere-tower v2: Q(ζ_12) + ∛ + √  (4·3·2 = 24).

Q(ζ_12): φ(12)=4, contains ζ_3 (rotation) → ∛ layer is Kummer C3.
Then √ layer on top. 圆⊕塔⊕塔.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260818
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40

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

def layer(f, bcoeffs, kth, do_polred=False):
    """Res_y(f(y), z^k − b(y)) → degree k·deg(f). bcoeffs ascending."""
    d = int(f.poldegree())
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs[::-1], Y)
    g = Z**kth - b
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    if R.poldegree() != kth * d:
        return None
    c = prim_part([int(z) for z in pari.Vecrev(pari.Pol(R, X))])
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    if do_polred and kth * d <= 12:
        try:
            Pr = pari.polredabs(P)
            if Pr.poldegree() == kth * d:
                P = Pr
        except Exception:
            pass
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return monic_ize(cc)

def rand_elem(f, rng, bound=4):
    d = int(f.poldegree())
    for _ in range(60):
        bc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def harm_elem(f, rng, n=12):
    """ζ_n-weighted sum (discrete harmonics)."""
    d = int(f.poldegree())
    for _ in range(60):
        w = [rng.randint(-2, 2) for _ in range(n)]
        if all(z == 0 for z in w):
            continue
        zeta = pari.Mod(X, f)
        acc = pari.Mod(0, f)
        zk = pari.Mod(1, f)
        for k in range(n):
            acc = acc + w[k] * zk
            zk = zk * zeta
        bc = [int(z) for z in pari.Vecrev(pari.Pol(pari.lift(acc), X))]
        while len(bc) < d:
            bc.append(0)
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def main():
    rng = random.Random(SEED)
    f0 = pari.polcyclo(12)  # Q(ζ12), degree 4, contains ζ3
    print(f"f0 = Q(ζ12) deg={f0.poldegree()}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 500 * N_TARGET:
        tries += 1
        bc = harm_elem(f0, rng) if rng.random() < 0.5 else rand_elem(f0, rng)
        if bc is None: continue
        c = layer(f0, bc, 3, do_polred=True)  # ∛ Kummer layer (ζ3 ∈ Q(ζ12))
        if c is None: continue
        f1 = pari.Polrev(c)  # degree 12
        bc = rand_elem(f1, rng, bound=3)
        if bc is None: continue
        c = layer(f1, bc, 2)  # √ layer
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
    print(f"sphere2: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # sphere2 r={r}")

if __name__ == "__main__":
    main()
