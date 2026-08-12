#!/usr/bin/env python3
"""IGP24 generator v7: degree-24 substitution families (wreath groups).

f(x) = g(x^k), deg g = 24/k, k in {2,3,4,6,8,12}:
  - k=2:  f = g(x^2), g deg 12 -> group C2 wr Gal(g)-ish; r(f)=2*#pos-roots(g)
  - k=3:  f = g(x^3), g deg 8;  r(f)=r(g)
  - k=4:  f = g(x^4), g deg 6;  r(f)=4*#pos-roots(g)
  - k=6:  f = g(x^6), g deg 4;  r(f)=6*#pos-roots(g) if k even
  - k=8:  f = g(x^8), g deg 3
  - k=12: f = g(x^12), g deg 2

Bases g are built diversely: perturbed-split (S_d / A_d / signature-controlled),
abelian subcyclo, wreath-of-smaller, product composita.
Target r buckets for degree-24 output: 24, 16, 20, 12, 8, 4, 0.
"""
import random, sys
from cypari2 import Pari

pari = Pari()
X = pari("x")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260818

def mk(c): return pari.Polrev(c)
def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1
def irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False
def coeffs(g):
    return [int(z) for z in pari.Vecrev(pari.Pol(g, X))]
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True
def mul_lin(poly, root):
    out = [0] * (len(poly) + 1)
    for i, ci in enumerate(poly):
        out[i] -= ci * root
        out[i + 1] += ci
    return out
def mul_quad(poly, p_):
    out = [0] * (len(poly) + 2)
    for i, ci in enumerate(poly):
        out[i] += p_ * ci
        out[i + 2] += ci
    return out

def split_base(d, n_real, rng, bound=10, eps_max=5, positive=False):
    if positive and n_real >= bound:
        bound = n_real + 8
    lo, hi = (1, bound) if positive else (-bound, bound)
    r_i = sorted(rng.sample(range(lo, hi + 1), n_real))
    if n_real >= 2 and min(r_i[i+1] - r_i[i] for i in range(n_real - 1)) < 1:
        return None
    n_pairs = (d - n_real) // 2
    p_j = [rng.randint(1, 8) for _ in range(n_pairs)]
    base = [1]
    for ri in r_i: base = mul_lin(base, ri)
    for pj in p_j: base = mul_quad(base, pj)
    for _ in range(30):
        k = rng.randint(0, d - 1)
        eps = rng.choice([z for z in range(-eps_max, eps_max + 1) if z != 0])
        c = base.copy(); c[k] += eps
        if c[0] == 0: continue
        g = mk(c)
        if not irred(g): continue
        if real_count(g) == n_real: return c
    return None

def subst(g, k):
    f = [0] * ((len(g) - 1) * k + 1)
    for i, gi in enumerate(g):
        f[k * i] = gi
    return f

def main():
    rng = random.Random(SEED)
    out = []
    seen = set()
    def emit(c, tag):
        key = tuple(c)
        if key in seen: return
        seen.add(key)
        out.append((c, tag))

    # k=2: deg-12 bases; want r24 = 2*#pos -> need 12 positive roots (totally real pos)
    #       r16 = 2*8 (8 pos roots), r12 = 2*6, r8 = 2*4, r4 = 2*2, r0 = 0 pos
    print("k=2 ...", file=sys.stderr)
    for npos, want_r in ((12, 24), (8, 16), (6, 12), (4, 8), (2, 4), (0, 0)):
        made = 0
        while made < 40:
            g = split_base(12, npos, rng, positive=(npos == 12))
            if g is None: continue
            f = subst(g, 2)
            if f[0] == 0: continue
            fg = mk(f)
            if not irred(fg): continue
            if real_count(fg) == want_r:
                emit(f, f"subst2 pos={npos}")
                made += 1
        print(f"  k=2 r={want_r}: {made}", file=sys.stderr)

    # k=3: deg-8 bases, r(f) = r(g); want r24: g totally real r=8... 8*3=24 no:
    #   r(f)=r(g) for odd k -> want r in {24} impossible (deg 8). r values: 0..8.
    print("k=3 ...", file=sys.stderr)
    for n_real, want_r in ((8, 8), (6, 6), (4, 4), (2, 2), (0, 0)):
        made = 0
        while made < 25:
            g = split_base(8, n_real, rng)
            if g is None: continue
            f = subst(g, 3)
            if f[0] == 0: continue
            fg = mk(f)
            if not irred(fg): continue
            if real_count(fg) == want_r:
                emit(f, f"subst3 r={n_real}")
                made += 1
        print(f"  k=3 r={want_r}: {made}", file=sys.stderr)

    # k=4: deg-6 bases; r(f) = 4*#pos(g) if all... x^4 = y has 2 real roots for y>0:
    #   r(f) = 2*(#positive real roots of g). want r24: 12 pos roots of deg-6 g:
    #   impossible -> max r = 2*6 = 12. r16 no. r12 = 6 pos, r8 = 4 pos, r4=2 pos, r0=0.
    print("k=4 ...", file=sys.stderr)
    for npos, want_r in ((6, 12), (4, 8), (2, 4), (0, 0)):
        made = 0
        while made < 25:
            g = split_base(6, npos, rng, positive=(npos == 6))
            if g is None: continue
            f = subst(g, 4)
            if f[0] == 0: continue
            fg = mk(f)
            if not irred(fg): continue
            if real_count(fg) == want_r:
                emit(f, f"subst4 pos={npos}")
                made += 1
        print(f"  k=4 r={want_r}: {made}", file=sys.stderr)

    # k=6: deg-4 bases; r(f) = 2*(#pos roots of g): r12 = 6 pos (deg4: max 4) no;
    #   r8 = 4 pos, r4 = 2 pos, r0 = 0.
    print("k=6 ...", file=sys.stderr)
    for npos, want_r in ((4, 8), (2, 4), (0, 0)):
        made = 0
        while made < 25:
            g = split_base(4, npos, rng, positive=(npos == 4))
            if g is None: continue
            f = subst(g, 6)
            if f[0] == 0: continue
            fg = mk(f)
            if not irred(fg): continue
            if real_count(fg) == want_r:
                emit(f, f"subst6 pos={npos}")
                made += 1
        print(f"  k=6 r={want_r}: {made}", file=sys.stderr)

    # k=8: deg-3 bases; r(f) = r(g) (odd k): r in {1,3}
    print("k=8 ...", file=sys.stderr)
    for n_real, want_r in ((3, 3), (1, 1)):
        made = 0
        while made < 20:
            g = split_base(3, n_real, rng)
            if g is None: continue
            f = subst(g, 8)
            if f[0] == 0: continue
            fg = mk(f)
            if not irred(fg): continue
            if real_count(fg) == want_r:
                emit(f, f"subst8 r={n_real}")
                made += 1
        print(f"  k=8 r={want_r}: {made}", file=sys.stderr)

    print(f"total: {len(out)}", file=sys.stderr)
    for c, tag in out:
        r = real_count(mk(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
