#!/usr/bin/env python3
"""升维塔 G: degree-12 totally-real base + ONE sqrt layer (12·2 = 24).

Block system of 24T1958 = 12 blocks × 2 points → field = K12(√b).
K12 = real cyclotomic subfield (φ(n)/2 = 12: n ∈ {35,39,45,52,56,70,72,78,84,90}).
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260815
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
N_BASE = int(sys.argv[3]) if len(sys.argv) > 3 else 35

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

def get_real_subfield(n, deg):
    res = pari.polsubcyclo(n, deg)
    polys = res if res.type() == "t_VEC" else [res]
    for p in polys:
        if p.poldegree() == deg and int(pari.polsturm(p)) == deg:
            try:
                pr = pari.polredabs(p)
                if pr.poldegree() == deg and int(pari.polsturm(pr)) == deg:
                    return pr
            except Exception:
                pass
            return p
    return None

def monic_ize(c):
    lc = c[-1]
    if lc == 1: return c
    if lc == -1: return [-z for z in c]
    if all(z % lc == 0 for z in c): return [z // lc for z in c]
    d = len(c) - 1
    out = [c[k] * lc ** (d - 1 - k) for k in range(d)]
    out.append(1)
    return out

def sqrt_layer(f, bcoeffs):
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
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return monic_ize(cc)

def b_lowdeg(f, rng, cmin, cmax, deg=2):
    d = int(f.poldegree())
    for _ in range(60):
        bc = [0] * d
        k = rng.randint(0, min(deg, d - 1))
        for j in range(1, k + 1):
            bc[j] = rng.randint(-1, 1)
        if all(z == 0 for z in bc):
            bc[0] = 1
        bc[0] += rng.randint(cmin, cmax)
        return bc
    return None

def main():
    rng = random.Random(SEED)
    f0 = get_real_subfield(N_BASE, 12)
    if f0 is None:
        print(f"FAIL: no deg-12 subfield of Q(ζ{N_BASE})", file=sys.stderr)
        return
    print(f"f0(Q(ζ{N_BASE})^+) deg={f0.poldegree()} r={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 400 * N_TARGET:
        tries += 1
        bc = b_lowdeg(f0, rng, 8, 12)  # positive shift → r=24-ish; use -4..8 for mix
        if bc is None:
            continue
        c = sqrt_layer(f0, bc)
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
    print(f"liftG: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # liftG n={N_BASE} r={r}")

if __name__ == "__main__":
    main()
