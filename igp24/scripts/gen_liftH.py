#!/usr/bin/env python3
"""升维塔 H: non-cyclotomic Abelian base (class field via quadray) + double √.

K0 = ray class field of Q(√D) mod m (Abelian, non-cyclotomic), degree 6.
L = K0(√b)(√c)  (6·2·2 = 24).

Kummer layers over an Abelian base are AUTOMATICALLY equivariant
(σ(ⁿ√a) = ζ^k·ⁿ√a by the image of ζ) — the "rotation" structure.  The
cyclotomic-real bases (ζ13 etc.) gave 24T21854; quadray bases give
different Abelian fields → different groups.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20261001
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
MODE = sys.argv[3] if len(sys.argv) > 3 else "full"  # full (r=24) | mix | neg
DISC = int(sys.argv[4]) if len(sys.argv) > 4 else -23
MODULUS = int(sys.argv[5]) if len(sys.argv) > 5 else 3

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

def quadray_base(D, m):
    """Ray class field over Q(√D) mod m; polredabs'd absolute poly."""
    try:
        res = pari.quadray(D, m)
        polys = res if res.type() == "t_VEC" else [res]
        best = None
        for p in polys:
            if p.type() != "t_POL" or p.poldegree() < 2:
                continue
            try:
                pr = pari.polredabs(p)
                if pr.poldegree() == p.poldegree():
                    return pr
            except Exception:
                pass
            return p
        return best
    except Exception:
        return None

def main():
    rng = random.Random(SEED)
    f0 = quadray_base(DISC, MODULUS)
    if f0 is None:
        print(f"FAIL: no quadray({DISC},{MODULUS}) base", file=sys.stderr)
        return
    print(f"f0 = quadray({DISC},{MODULUS}) deg={f0.poldegree()} r={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 400 * N_TARGET:
        tries += 1
        if MODE == "full":
            bc = b_lowdeg(f0, rng, 8, 12)
        elif MODE == "neg":
            bc = b_lowdeg(f0, rng, -12, -8)
        else:
            bc = b_lowdeg(f0, rng, -4, 8)
        if bc is None: continue
        c1 = sqrt_layer(f0, bc, do_polred=True)
        if c1 is None: continue
        f1 = pari.Polrev(c1)
        if MODE == "full":
            cc = b_lowdeg(f1, rng, 8, 12)
        elif MODE == "neg":
            cc = b_lowdeg(f1, rng, -12, -8)
        else:
            cc = b_lowdeg(f1, rng, -4, 8)
        if cc is None: continue
        c2 = sqrt_layer(f1, cc)
        if c2 is None: continue
        if len(c2) != 25: continue
        if c2[0] == 0: continue
        key = tuple(c2)
        if key in seen: continue
        g = pari.Polrev(c2)
        if not is_irred(g): continue
        r = int(pari.polsturm(g))
        seen.add(key)
        out.append((c2, r))
        made += 1
    print(f"liftH_{MODE}: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # liftH_{MODE} r={r}")

if __name__ == "__main__":
    main()
