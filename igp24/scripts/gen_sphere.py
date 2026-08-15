#!/usr/bin/env python3
"""Sphere-tower: cyclotomic (Euler φ) base + Kummer tower layer (圆⊕塔).

L = Q(ζ_15)(∛a), a ∈ Q(ζ_15):  [L:Q] = 8·3 = 24.
- 圆层: Q(ζ_15), φ(15) = 8 (Gal (Z/15)* ≅ C2×C4 — Euler rotation group)
- 塔层: ∛a over Q(ζ_15) — Kummer C3 (ζ_3 ∈ Q(ζ_15) since 3 | 15)
- ζ_3 ∈ Q(ζ_15) makes the cubic layer cyclic (rotation by ζ_3^k)
- Sphere link: a chosen as symmetric sums of ζ_15 powers (discrete harmonics)
  or generic elements; both variants generated.

Variants:
  BASE=15  (cyclotomic, 8·3)
  BASE=21  (cyclotomic, φ(21)=12, ∛ over degree-12? no: 12·3=36 ≠ 24 — skip)
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260817
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
MODE = sys.argv[3] if len(sys.argv) > 3 else "gen"  # gen | harm (harmonic symmetric sums)
N_BASE = int(sys.argv[4]) if len(sys.argv) > 4 else 15  # cyclotomic order (φ(n)=8: 15,16,20,24,30)

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

def cubic_layer(f, bcoeffs):
    """Res_y(f(y), z³ − b(y)) → degree 3·deg(f). bcoeffs ascending."""
    d = int(f.poldegree())
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs[::-1], Y)
    g = Z**3 - b
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    if R.poldegree() != 3 * d:
        return None
    c = prim_part([int(z) for z in pari.Vecrev(pari.Pol(R, X))])
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
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

def harm_elem(f, rng, n=15):
    """Sphere/harmonic element: b = Σ ζ_n^{k}·c_k — symmetric sums of roots of unity
    (discrete spherical harmonics on the cyclotomic lattice)."""
    d = int(f.poldegree())
    # b = ζ15 + ζ15^2·c1 + ... : use coefficients as "weights" on roots of unity
    for _ in range(60):
        w = [rng.randint(-2, 2) for _ in range(n)]
        if all(z == 0 for z in w):
            continue
        # sum w_k · ζ15^k mod f(ζ15)
        # represent ζ15^k in the power basis via Mod
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
    f0 = pari.polcyclo(N_BASE)  # Q(ζ_N_BASE), degree φ(n) = 8
    print(f"f0 = Q(ζ{N_BASE}) deg={f0.poldegree()}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 400 * N_TARGET:
        tries += 1
        if MODE == "harm":
            bc = harm_elem(f0, rng)
        else:
            bc = rand_elem(f0, rng)
        if bc is None:
            continue
        c = cubic_layer(f0, bc)
        if c is None:
            continue
        if len(c) != 25:
            continue
        if c[0] == 0:
            continue
        key = tuple(c)
        if key in seen:
            continue
        g = pari.Polrev(c)
        if not is_irred(g):
            continue
        r = int(pari.polsturm(g))
        seen.add(key)
        out.append((c, r))
        made += 1
    print(f"sphere({MODE}): made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # sphere_{MODE} r={r}")

if __name__ == "__main__":
    main()
