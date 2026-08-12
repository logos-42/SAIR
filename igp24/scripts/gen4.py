#!/usr/bin/env python3
"""IGP24 generator v4: diverse small-degree base pools -> composita.

Base-pool diversity (degrees 3,4,6,8,12):
  - abelian: polsubcyclo(p, d) real subfields (C3,C4,C6,C8,C12), x^d-a (C_d)
  - dihedral-ish: x^6+ax^3+b (D6), biquadratics (V4/D4/C4)
  - wreath: g(x^k) for deg-g * k = d
  - product: composita of smaller degrees
  - generic: perturbed split (S_d), random (S_d)
Composita over (a,b) with a*b=24, targeting open r buckets.
"""
import random, sys
from cypari2 import Pari

pari = Pari()
X = pari("x")

def mk(c):
    return pari.Polrev(c)

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

def split_base(d, n_real, rng, bound=12, eps_max=6):
    r_i = sorted(rng.sample(range(-bound, bound + 1), n_real))
    if n_real >= 2 and min(r_i[i+1] - r_i[i] for i in range(n_real - 1)) < 1:
        return None
    n_pairs = (d - n_real) // 2
    p_j = [rng.randint(1, 9) for _ in range(n_pairs)]
    base = [1]
    for ri in r_i: base = mul_lin(base, ri)
    for pj in p_j: base = mul_quad(base, pj)
    for _ in range(25):
        k = rng.randint(0, d - 1)
        eps = rng.choice([z for z in range(-eps_max, eps_max + 1) if z != 0])
        c = base.copy(); c[k] += eps
        if c[0] == 0: continue
        g = mk(c)
        if not irred(g): continue
        if real_count(g) == n_real: return c
    return None

def subcyclo(p, d):
    """real subfield of Q(zeta_p) of degree d (totally real, cyclic)."""
    try:
        res = pari.polsubcyclo(p, d)
    except Exception:
        return None
    items = list(res) if res.type() == "t_VEC" else [res]
    for f in items:
        try:
            if f.poldegree() == d and irred(f) and real_count(f) == d:
                return coeffs(f)
        except Exception:
            continue
    return None

def subst(g, k):
    f = [0] * (len(g) - 1) * k + [0] * k  # placeholder
    f = [0] * ((len(g) - 1) * k + 1)
    for i, gi in enumerate(g):
        f[k * i] = gi
    return f

def compositum(fc, gc):
    try:
        comps = pari.polcompositum(mk(fc), mk(gc))
    except Exception:
        return []
    deg = (len(fc) - 1) * (len(gc) - 1)
    out = []
    for R in comps:
        try:
            if R.poldegree() == deg and irred(R):
                out.append(coeffs(R))
        except Exception:
            continue
    return out

def build_pool(d, rng, n_each=30):
    """Build pool of degree-d polys with varied groups; return list of (coeffs, n_real)."""
    pool = []
    seen = set()

    def add(c, maxn=60):
        if c is None or len(c) - 1 != d: return
        key = tuple(c)
        if key in seen or len(pool) >= maxn: return
        g = mk(c)
        if not irred(g): return
        r = real_count(g)
        if r < 0: return
        seen.add(key)
        pool.append((c, r))

    if d == 3:
        # S3 cubic (totally real disc>0, r=1 disc<0), A3
        for a in range(1, 10):
            for b in range(1, 10):
                c = [-b, -a, 0, 1]
                add(c)
                if len(pool) >= 40: break
            if len(pool) >= 40: break
        for c in ([1, -3, 0, 1], [-1, 3, 0, 1]):
            add(c)
    elif d == 4:
        # biquadratics x^4 + a x^2 + b (D4/V4/C4), x^4 - a (C4)
        for a in range(-14, 15):
            if a == 0: continue
            for b in range(1, 12):
                c = [b, 0, a, 0, 1]
                add(c)
        for a in range(2, 14):
            c = [-a, 0, 0, 0, 1]
            add(c)
        # generic S4/A4
        for _ in range(60):
            c = split_base(4, rng.choice((2, 4)), rng)
            add(c)
    elif d == 6:
        # x^6 + a x^3 + b (D6), x^6 - a (C6)
        for a in range(-8, 9):
            for b in range(1, 9):
                c = [b, 0, 0, a, 0, 0, 1]
                add(c)
        for a in range(2, 14):
            c = [-a, 0, 0, 0, 0, 0, 1]
            add(c)
        # polsubcyclo C6 real
        p = 31
        while p < 300 and len(pool) < 10:
            if is_prime(p) and (p - 1) % 12 == 0:
                add(subcyclo(p, 6))
            p += 12
        # wreath g(x^2), g deg 3
        for _ in range(40):
            g = split_base(3, rng.choice((1, 3)), rng)
            if g: add(subst(g, 2))
        # composita (2,3)
        for _ in range(40):
            q = split_base(2, 2, rng)
            g = split_base(3, rng.choice((1, 3)), rng)
            if q and g:
                for c in compositum(q, g)[:2]:
                    add(c)
    elif d == 8:
        # polsubcyclo C8 real
        p = 17
        while p < 500 and len(pool) < 10:
            if is_prime(p) and (p - 1) % 16 == 0:
                add(subcyclo(p, 8))
            p += 16
        # wreath g(x^2) deg 4, g(x^4) deg 2
        for _ in range(50):
            g = split_base(4, rng.choice((0, 2, 4)), rng)
            if g: add(subst(g, 2))
        for _ in range(30):
            g = split_base(2, 2, rng)
            if g: add(subst(g, 4))
        # composita (4,2),(2,4)
        for _ in range(50):
            q = split_base(2, 2, rng)
            g = split_base(4, rng.choice((0, 2, 4)), rng)
            if q and g:
                for c in compositum(q, g)[:2]:
                    add(c)
        # generic
        for _ in range(40):
            c = split_base(8, rng.choice((0, 2, 4, 6, 8)), rng)
            add(c)
    elif d == 12:
        # polsubcyclo C12 real
        p = 37
        while p < 600 and len(pool) < 12:
            if is_prime(p) and (p - 1) % 24 == 0:
                add(subcyclo(p, 12))
            p += 24
        # wreath: g(x^2) deg 6, g(x^3) deg 4, g(x^4) deg 3, g(x^6) deg 2
        for k, gd in ((2, 6), (3, 4), (4, 3), (6, 2)):
            for _ in range(25):
                g = split_base(gd, rng.choice((0, gd if gd % 2 == 0 else 1)), rng)
                if g: add(subst(g, k))
        # composita (6,2),(4,3),(3,4),(2,6)
        for (a, b) in ((6, 2), (4, 3), (3, 4), (2, 6)):
            for _ in range(20):
                pa = split_base(a, rng.choice((a if a % 2 == 0 else 1, 0, 2)), rng)
                pb = split_base(b, rng.choice((b if b % 2 == 0 else 1, 0, 2)), rng)
                if pa and pb:
                    for c in compositum(pa, pb)[:2]:
                        add(c)
        # generic totally real
        for _ in range(30):
            c = split_base(12, rng.choice((0, 4, 8, 12)), rng)
            add(c)
    return pool

def main():
    rng = random.Random(20260816)
    pools = {}
    for d in (3, 4, 6, 8, 12):
        pools[d] = build_pool(d, rng)
        print(f"pool d={d}: {len(pools[d])}  (r-dist: "
              f"{ {r: sum(1 for _, rr in pools[d] if rr == r) for r in sorted(set(rr for _, rr in pools[d]))} })",
              file=sys.stderr)

    out = []
    seen = set()

    def emit(c, tag):
        key = tuple(c)
        if key in seen: return
        seen.add(key)
        out.append((c, tag))

    pairs = [(2,12),(3,8),(4,6),(6,4),(8,3),(12,2)]
    for want_r, n in ((24, 320), (16, 200), (20, 150), (12, 100), (8, 80),
                      (4, 60), (0, 60), (18, 40), (6, 40)):
        made = 0; tries = 0
        while made < n and tries < 5000:
            tries += 1
            a, b = rng.choice(pairs)
            pa = pools.get(a); pb = pools.get(b)
            if not pa or not pb: continue
            fa, r1 = rng.choice(pa)
            gb, r2 = rng.choice(pb)
            if r1 * r2 != want_r: continue
            fs = compositum(fa, gb)
            if not fs: continue
            for cf in fs:
                if real_count(mk(cf)) != want_r: continue
                emit(cf, f"comp {a}x{b} ({r1}*{r2})")
                made += 1
                break
        print(f"composita r={want_r}: {made}", file=sys.stderr)

    print(f"total: {len(out)}", file=sys.stderr)
    for c, tag in out:
        r = real_count(mk(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
