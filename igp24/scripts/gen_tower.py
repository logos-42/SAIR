#!/usr/bin/env python3
"""Tower generator v0: K8(∛b) tower constructions for solvable 2-3 groups.

L = K8(∛b), K8 = Q(√a1,√a2,√a3) triquadratic (degree 8), b ∈ K8 random.
Minimal poly of ∛b over Q = Res_y(f8(y), z³−b(y)) with z→x: degree 24.
Random b sweeps solvable 2-3 groups (order-768 family etc).
"""
import math, random, sys
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260819

def squarefree_int(rng, bound=60, allow_neg=True):
    while True:
        a = rng.randint(2, bound)
        if a % 4 == 0: continue
        sq = True
        d = 2
        while d * d <= a:
            if a % (d * d) == 0:
                sq = False
                break
            d += 1
        if sq:
            return -a if (allow_neg and rng.random() < 0.5) else a

def is_irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False

def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1

def triquadratic(a1, a2, a3):
    f2 = pari.Polrev([-a1, 0, 1])
    c = pari.polcompositum(f2, pari.Polrev([-a2, 0, 1]))
    f4 = c[0]
    if f4.poldegree() != 4: return None
    c = pari.polcompositum(f4, pari.Polrev([-a3, 0, 1]))
    f8 = c[0]
    if f8.poldegree() != 8: return None
    if not is_irred(f8): return None
    return f8

def tower_poly(f8, bc):
    """Res_y(f8(y), z³−b(y)) with z→x; primitive part of the degree-24 poly."""
    f8y = pari.subst(f8, X, Y)
    coeffs = [Z**3 - bc[0]] + [-bc[i] for i in range(1, len(bc))]
    g = pari.Pol(coeffs, Y)
    try:
        R = pari.polresultant(f8y, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    if R.poldegree() != 24:
        return None
    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
    # primitive part
    g = 0
    for z in c:
        g = math.gcd(g, abs(z))
    if g > 1:
        c = [z // g for z in c]
    return c

def monicize(c):
    """Rule formula: g_k = a_k * a_24^(23-k)  (monic, same number field)."""
    lc = c[24]
    if lc == 1:
        return c
    if lc == -1:
        return [-z for z in c]
    if all(z % lc == 0 for z in c):
        return [z // lc for z in c]
    out = []
    for k in range(25):
        out.append(c[k] * lc ** (23 - k))
    return out

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    while made < n_target and tries < 6000:
        tries += 1
        a1 = squarefree_int(rng)
        a2 = squarefree_int(rng)
        a3 = squarefree_int(rng)
        if len({a1, a2, a3}) < 3: continue
        f8 = triquadratic(a1, a2, a3)
        if f8 is None: continue
        for _ in range(6):
            bc = [rng.randint(-4, 4) for _ in range(8)]
            if all(c == 0 for c in bc): continue
            c = tower_poly(f8, bc)
            if c is None: continue
            c = monicize(c)
            key = tuple(c)
            if key in seen: continue
            if c[0] == 0: continue
            g = pari.Polrev(c)
            if not is_irred(g): continue
            r = real_count(g)
            if r < 0 or r % 2 != 0: continue
            seen.add(key)
            out.append((c, f"tower a=({a1},{a2},{a3})"))
            made += 1
            break
    print(f"tower: made={made} tries={tries}", file=sys.stderr)
    for c, tag in out:
        r = real_count(pari.Polrev(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
