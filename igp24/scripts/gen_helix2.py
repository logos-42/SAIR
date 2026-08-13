#!/usr/bin/env python3
"""Helix v2 (totally-real): cyclotomic-style cubic + positive Kummer layers.

L0 = totally real cubic (disc>0, S3 or C3 Gal)
L1 = L0(√(γ²+1)),  γ ∈ L0   (all conjugates ≥ 1 → stays totally real)
L2 = L1(√(γ²+1)),  γ ∈ L1
L3 = L2(√(γ²+1)),  γ ∈ L2
degree 3·2·2·2 = 24, r(L) = 24.

Rotation (conjugation) acts on γ; positive elements keep the helix real.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260823

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

def pos_layer(f, gamma_coeffs, c, do_reduce=True):
    """Res_y(f(y), z²−(γ(y)²+c)) with z→x; monic degree-2d poly. c>0 rational."""
    d = f.poldegree()
    fy = pari.subst(f, X, Y)
    gam = pari.Pol(gamma_coeffs, Y)
    b = (gam * gam + c) % fy
    bc = [int(z) for z in pari.Vecrev(pari.Pol(b, Y))]
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
    c2 = prim_part([int(z) for z in pari.Vecrev(pari.Pol(R, X))])
    P = pari.Polrev(c2)
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

def rand_gamma(f, rng, bound=2):
    d = f.poldegree()
    for _ in range(60):
        gc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in gc):
            continue
        return gc
    return None

def totreal_cubic(rng, bound=7):
    for _ in range(400):
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
    n_target = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    while made < n_target and tries < 1500:
        tries += 1
        c3 = totreal_cubic(rng)
        if c3 is None: continue
        f = pari.Polrev(c3)
        gc = rand_gamma(f, rng)
        if gc is None: continue
        c = pos_layer(f, gc, rng.choice((1, 2, 3, 5)))
        if c is None: continue
        f2 = pari.Polrev(c)
        gc = rand_gamma(f2, rng)
        if gc is None: continue
        c = pos_layer(f2, gc, rng.choice((1, 2, 3, 5)))
        if c is None: continue
        f3 = pari.Polrev(c)
        gc = rand_gamma(f3, rng)
        if gc is None: continue
        c = pos_layer(f3, gc, rng.choice((1, 2, 3, 5)), do_reduce=False)
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
        if r != 24: continue
        seen.add(key)
        out.append((c, f"helix2 cube={tuple(c3)}"))
        made += 1
    print(f"helix2: made={made} tries={tries}", file=sys.stderr)
    for c, tag in out:
        r = real_count(pari.Polrev(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
