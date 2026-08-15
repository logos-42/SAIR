#!/usr/bin/env python3
"""升维塔 F: degree-8 totally-real base + cubic layer x³+px+q (disc>0).

K0 = Q(ζ17)^+ (degree 8, totally real, Gal C8)
L = K0(α), α root of x³ + p x + q, p,q ∈ K0, disc = −4p³−27q² > 0 (totally real cubic)
→ [L:Q] = 24, r(L) = 3·8 = 24 (each real embedding of K0 gets 3 real roots).

The cubic layer is NOT Kummer (ζ3 ∉ K0) → S3-type Galois data → different
group family than the √-towers (24T21854/22794). Targets r=24 open pairs.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260815
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
N_BASE = int(sys.argv[3]) if len(sys.argv) > 3 else 17  # cyclotomic order for degree-8 real subfield

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

def cubic_layer(f, pcoeffs, qcoeffs):
    """Res_y(f(y), z³ + p(y)·z + q(y)) → degree 3·deg(f) poly."""
    d = int(f.poldegree())
    fy = pari.subst(f, X, Y)
    # g as poly in y, coeffs in Z[z]:
    #   g = p(y)·z + q(y) + z³
    #   y^0: q0 + z³ ; y^i: qi + pi·z
    coeffs = [Z**3 + qcoeffs[0]] + [qcoeffs[i] + pcoeffs[i] * Z for i in range(1, d)]
    g = pari.Pol(coeffs, Y)
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
    lc = cc[-1]
    if lc == 1:
        return cc
    if lc == -1:
        return [-z for z in cc]
    if all(z % lc == 0 for z in cc):
        return [z // lc for z in cc]
    dd = len(cc) - 1
    out = [cc[k] * lc ** (dd - 1 - k) for k in range(dd)]
    out.append(1)
    return out

def rand_pq(f, rng, neg_shift=(-6, -3)):
    """p = −(γ²+c) (all-negative conjugates → −4p³ > 0), q small random."""
    d = int(f.poldegree())
    for _ in range(80):
        gc = [rng.randint(-1, 1) for _ in range(min(3, d))] + [0] * (d - min(3, d))
        if all(z == 0 for z in gc):
            continue
        gam = pari.Polrev(gc, X)
        c = rng.choice((1, 2, 3))
        p = -(gam * gam + c) % f
        pc = [int(z) for z in pari.Vecrev(pari.Pol(p, X))]
        while len(pc) < d:
            pc.append(0)
        qc = [rng.randint(-2, 2) for _ in range(d)]
        qc[0] = rng.randint(-2, 2)
        return pc, qc
    return None, None

def main():
    rng = random.Random(SEED)
    f0 = get_real_subfield(N_BASE, 8)
    if f0 is None:
        print(f"FAIL: no Q(ζ{N_BASE})^+ base", file=sys.stderr)
        return
    print(f"f0(Q(ζ{N_BASE})^+) deg={f0.poldegree()} r={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 500 * N_TARGET:
        tries += 1
        pc, qc = rand_pq(f0, rng)
        if pc is None:
            continue
        c = cubic_layer(f0, pc, qc)
        if c is None:
            continue
        key = tuple(c)
        if key in seen:
            continue
        if c[0] == 0:
            continue
        g = pari.Polrev(c)
        if not is_irred(g):
            continue
        r = real_count(g)
        if r < 0 or r % 2 != 0:
            continue
        seen.add(key)
        out.append((c, r))
        made += 1
    print(f"liftF: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # liftF r={r}")

if __name__ == "__main__":
    main()
