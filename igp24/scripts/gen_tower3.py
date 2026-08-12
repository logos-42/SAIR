#!/usr/bin/env python3
"""Tower generator v3: totally-real nested quadratic tower (r=24 target).

L = K3(√b1)(√b2)(√b3), K3 = totally real cubic (disc>0), b_i = γ_i² + 1 ∈ K
(all conjugates ≥ 1 → every layer keeps all roots real → r(L) = 24).
Nested √-tower with γ random gives varied non-abelian 2-groups ⋊ (S3/A3 data).

Degree: 3·2·2·2 = 24. Output: monic degree-24 polys with r=24.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260821

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

def sq_layer(f, gamma_coeffs, do_reduce=True):
    """Res_y(f(y), z²−(γ(y)²+1)) with z→x; monic degree-2d poly."""
    d = f.poldegree()
    fy = pari.subst(f, X, Y)
    gam = pari.Pol(gamma_coeffs, Y)
    b = (gam * gam + 1) % fy   # γ²+1 mod f
    bc = [int(z) for z in pari.Vecrev(pari.Pol(b, Y))]
    # pad to d coeffs
    while len(bc) < d:
        bc.append(0)
    coeffs = [Z**2 - bc[0]] + [-bc[i] for i in range(1, d)]
    g = pari.Pol(coeffs, Y)
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception:
        return None
    R = pari.subst(R, Z, X)
    d2 = 2 * d
    if R.poldegree() != d2:
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

def random_gamma(f, rng, bound=4):
    d = f.poldegree()
    for _ in range(30):
        gc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in gc):
            continue
        return gc
    return None

def totreal_cubic(rng, bound=8):
    """Random totally real irreducible cubic (disc > 0)."""
    for _ in range(300):
        c = [rng.randint(-bound, bound), rng.randint(-bound, bound), 0, 1]
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        if int(pari.poldisc(g)) > 0 and real_count(g) == 3:
            return c
    return None

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    made = 0
    tries = 0
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 150
    while made < n_target and tries < 2000:
        tries += 1
        c3 = totreal_cubic(rng)
        if c3 is None: continue
        f = pari.Polrev(c3)
        gc = random_gamma(f, rng)
        if gc is None: continue
        c = sq_layer(f, gc)
        if c is None: continue
        f2 = pari.Polrev(c)
        gc = random_gamma(f2, rng, bound=3)
        if gc is None: continue
        c = sq_layer(f2, gc)
        if c is None: continue
        f3 = pari.Polrev(c)
        gc = random_gamma(f3, rng, bound=3)
        if gc is None: continue
        c = sq_layer(f3, gc, do_reduce=False)
        if c is None: continue
        if len(c) != 25: continue
        c = monicize(c)
        key = tuple(c)
        if key in seen: continue
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        r = real_count(g)
        if r != 24:
            continue
        seen.add(key)
        out.append((c, f"tower3 cube={tuple(c3)}"))
        made += 1
    print(f"tower3: made={made} tries={tries}", file=sys.stderr)
    for c, tag in out:
        r = real_count(pari.Polrev(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
