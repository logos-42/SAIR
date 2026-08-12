#!/usr/bin/env python3
"""IGP24 generator v5: diverse bases with correct r coverage -> composita.

Fixes vs gen4: d=2 pool, totally-real-first pool ordering, positive-root bases
for x^2-substitution, higher pool caps.
"""
import random, sys
from cypari2 import Pari

pari = Pari()
X = pari("x")

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20260817

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

def split_base(d, n_real, rng, bound=12, eps_max=6, positive=False):
    """Irreducible deg-d poly with n_real real roots (or all-positive roots)."""
    lo, hi = (1, bound) if positive else (-bound, bound)
    r_i = sorted(rng.sample(range(lo, hi + 1), n_real))
    if n_real >= 2 and min(r_i[i+1] - r_i[i] for i in range(n_real - 1)) < 1:
        return None
    n_pairs = (d - n_real) // 2
    p_j = [rng.randint(1, 9) for _ in range(n_pairs)]
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

def subcyclo(p, d):
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

class Pool:
    def __init__(self, d, cap=200):
        self.d = d; self.cap = cap
        self.items = []  # (coeffs, n_real)
        self.seen = set()
    def add(self, c):
        if c is None or len(c) - 1 != self.d: return
        key = tuple(c)
        if key in self.seen or len(self.items) >= self.cap: return
        g = mk(c)
        if not irred(g): return
        r = real_count(g)
        if r < 0: return
        self.seen.add(key)
        self.items.append((c, r))

def build_pools(rng):
    pools = {d: Pool(d) for d in (2, 3, 4, 6, 8, 12)}
    # ---- degree 2: x^2 - a ----
    for a in range(2, 80):
        if a % 4 == 0: continue
        c = [-a, 0, 1]
        pools[2].add(c)
        c = [a, 0, 1]
        pools[2].add(c)
    # ---- degree 3 ----
    for a in range(1, 12):
        for b in range(1, 12):
            pools[3].add([-b, -a, 0, 1])   # S3, r=3 or 1
    for c in ([1, -3, 0, 1], [-1, 3, 0, 1]):  # A3
        pools[3].add(c)
    for _ in range(40):
        pools[3].add(split_base(3, rng.choice((1, 3)), rng))
    # ---- degree 4 ----
    for _ in range(30):
        pools[4].add(split_base(4, rng.choice((0, 2, 4)), rng))
    for a in range(-14, 15):
        if a == 0: continue
        for b in range(1, 12):
            pools[4].add([b, 0, a, 0, 1])   # biquadratic r=4 (a<0,a^2>4b) or r=0
    for a in range(2, 14):
        pools[4].add([-a, 0, 0, 0, 1])      # x^4-a, C4, r=2 (a>0) / r=0
    # ---- degree 6 ----
    for _ in range(35):
        pools[6].add(split_base(6, rng.choice((0, 2, 4, 6)), rng))
    for a in range(-8, 9):
        if a == 0: continue
        for b in range(1, 9):
            pools[6].add([b, 0, 0, a, 0, 0, 1])   # x^6+ax^3+b: r=6 if a<0,b>0,a^2>4b
    for a in range(2, 14):
        pools[6].add([-a, 0, 0, 0, 0, 0, 1])      # C6
    p = 31
    while p < 300:
        if is_prime(p) and (p - 1) % 12 == 0:
            pools[6].add(subcyclo(p, 6))
        p += 12
    for _ in range(30):
        g = split_base(3, rng.choice((1, 3)), rng)
        if g: pools[6].add(subst(g, 2))
    for _ in range(30):
        q = split_base(2, 2, rng); g = split_base(3, rng.choice((1, 3)), rng)
        if q and g:
            for c in compositum(q, g)[:2]: pools[6].add(c)
    # ---- degree 8 ----
    for _ in range(35):
        pools[8].add(split_base(8, rng.choice((0, 2, 4, 6, 8)), rng))
    p = 17
    while p < 500:
        if is_prime(p) and (p - 1) % 16 == 0:
            pools[8].add(subcyclo(p, 8))
        p += 16
    for _ in range(30):
        g = split_base(4, rng.choice((0, 2, 4)), rng)
        if g: pools[8].add(subst(g, 2))
    for _ in range(20):
        q = split_base(2, 2, rng); g = split_base(4, rng.choice((0, 2, 4)), rng)
        if q and g:
            for c in compositum(q, g)[:2]: pools[8].add(c)
    # ---- degree 12 ----
    for _ in range(35):
        pools[12].add(split_base(12, rng.choice((0, 2, 4, 6, 8, 10, 12)), rng))
    for _ in range(40):
        pools[12].add(split_base(12, 10, rng))   # r=10 bases for r=20 composita
    p = 37
    while p < 600:
        if is_prime(p) and (p - 1) % 24 == 0:
            pools[12].add(subcyclo(p, 12))
        p += 24
    for k, gd in ((2, 6), (3, 4), (4, 3), (6, 2)):
        for _ in range(20):
            g = split_base(gd, rng.choice((0, gd if gd % 2 == 0 else 1, 2)), rng)
            if g: pools[12].add(subst(g, k))
    for (a, b) in ((6, 2), (4, 3), (3, 4), (2, 6)):
        for _ in range(15):
            pa = split_base(a, rng.choice((a if a % 2 == 0 else 1, 0, 2)), rng)
            pb = split_base(b, rng.choice((b if b % 2 == 0 else 1, 0, 2)), rng)
            if pa and pb:
                for c in compositum(pa, pb)[:2]: pools[12].add(c)
    return pools

def main():
    rng = random.Random(SEED)
    pools = build_pools(rng)
    for d, p in pools.items():
        dist = {}
        for _, r in p.items:
            dist[r] = dist.get(r, 0) + 1
        print(f"pool d={d}: {len(p.items)} {dist}", file=sys.stderr)

    out = []
    seen = set()
    def emit(c, tag):
        key = tuple(c)
        if key in seen: return
        seen.add(key)
        out.append((c, tag))

    pairs = [(2,12),(3,8),(4,6),(6,4),(8,3),(12,2)]
    for want_r, n in ((24, 300), (16, 200), (20, 150), (12, 100), (8, 80),
                      (4, 60), (0, 60), (18, 40), (6, 40)):
        made = 0; tries = 0
        while made < n and tries < 6000:
            tries += 1
            a, b = rng.choice(pairs)
            pa = pools[a]; pb = pools[b]
            if not pa.items or not pb.items: continue
            fa, r1 = rng.choice(pa.items)
            gb, r2 = rng.choice(pb.items)
            if r1 * r2 != want_r: continue
            fs = compositum(fa, gb)
            if not fs: continue
            for cf in fs:
                if real_count(mk(cf)) != want_r: continue
                emit(cf, f"comp {a}x{b} ({r1}*{r2})")
                made += 1
                break
        print(f"composita r={want_r}: {made}", file=sys.stderr)

    # direct degree-24: x^2-substitution of positive-root deg-12 bases (totally real, wreath groups)
    print("x^2-substitution...", file=sys.stderr)
    pos12 = [c for c, r in pools[12].items if r == 12]
    if not pos12:
        for _ in range(60):
            c = split_base(12, 12, rng, positive=True)
            if c: pos12.append(c)
    for c in pos12:
        f = subst(c, 2)
        if f[0] == 0: continue
        g = mk(f)
        if irred(g) and real_count(g) == 24:
            emit(f, "subst2 pos12")

    print(f"total: {len(out)}", file=sys.stderr)
    for c, tag in out:
        r = real_count(mk(c))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
