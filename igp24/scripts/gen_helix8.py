#!/usr/bin/env python3
"""Helix v8 (K3-element sqrt tower): all sqrt elements live in K3.

K3 = totally-real S3 cubic.  L = K3(√a)(√b)(√c) with a, b, c ∈ K3
(NOT in the upper layers).  Conjugates of each layer element are governed
by the S3 Galois orbit of K3 → the closure is small (target: 24T1958 =
(D8:(C4×C4)):S3, order 768 = 6·128; S3 layer from K3, 2-group from the
correlated √ elements).

Difference vs helix5: helix5 drew √ elements randomly from the CURRENT
layer (independent conjugates → closure 2^24); here they come from K3
(S3-correlated → small 2-part).  r control: choose a,b,c totally positive
(positive shift) for r=24, or mixed for diversity.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X, Y, Z = pari("x"), pari("y"), pari("z")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260915
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
MODE = sys.argv[3] if len(sys.argv) > 3 else "mix"  # full (r=24) | mix | neg

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

def k3_elem(rng, mode, shift=(6, 10)):
    """Element of K3 (coefficients in the power basis): positive shift for
    totally-positive (r=24), negative for r=0, mixed otherwise."""
    for _ in range(60):
        if mode == "full":
            bc = [rng.randint(shift[0], shift[1])] + [rng.randint(-1, 1) for _ in range(2)]
        elif mode == "neg":
            bc = [-rng.randint(shift[0], shift[1])] + [rng.randint(-1, 1) for _ in range(2)]
        else:
            bc = [rng.randint(-shift[1], shift[1])] + [rng.randint(-2, 2) for _ in range(2)]
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def embed_k3_elem(bc, f_base, f_upper):
    """Embed element b(θ) ∈ K3 into K1 (upper field) via nfisincl.

    Returns coefficient vector of b in the basis of the upper field."""
    try:
        imgs = pari.nfisincl(f_base, f_upper)
    except Exception:
        return None
    if not imgs:
        return None
    img = imgs[0]  # θ ↦ img(β) polynomial in the upper field's generator
    bpoly = pari.Polrev(bc, X)
    b_up = pari.subst(bpoly, X, img) % f_upper
    cu = [int(z) for z in pari.Vecrev(pari.Pol(b_up, X))]
    d = int(f_upper.poldegree())
    while len(cu) < d:
        cu.append(0)
    return cu

def s3_cubic(rng, bound=7):
    for _ in range(400):
        c = [rng.randint(-bound, bound), rng.randint(-bound, bound), 0, 1]
        if c[0] == 0: continue
        g = pari.Polrev(c)
        if not is_irred(g): continue
        disc = int(pari.poldisc(g))
        if disc > 0 and real_count(g) == 3 and not pari.issquare(disc):
            return c
    return None

def main():
    rng = random.Random(SEED)
    c3 = s3_cubic(rng)
    if c3 is None:
        print("FAIL: no S3 cubic", file=sys.stderr)
        return
    f0 = pari.Polrev(c3)
    print(f"K3: {c3} disc={int(pari.poldisc(f0))} r={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 400 * N_TARGET:
        tries += 1
        # three layers, ALL elements from K3 (embedded via nfisincl)
        bc = k3_elem(rng, MODE)
        if bc is None: continue
        c = sqrt_layer(f0, bc, do_polred=True)
        if c is None: continue
        f1 = pari.Polrev(c)  # degree 6
        bc = k3_elem(rng, MODE)
        if bc is None: continue
        cu = embed_k3_elem(bc, f0, f1)
        if cu is None: continue
        c = sqrt_layer(f1, cu, do_polred=True)
        if c is None: continue
        f2 = pari.Polrev(c)  # degree 12
        bc = k3_elem(rng, MODE)
        if bc is None: continue
        cu = embed_k3_elem(bc, f0, f2)
        if cu is None: continue
        c = sqrt_layer(f2, cu)
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
    print(f"helix8_{MODE}: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # helix8_{MODE} r={r}")

if __name__ == "__main__":
    main()
