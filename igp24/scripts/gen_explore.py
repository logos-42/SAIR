#!/usr/bin/env python3
"""IGP24 candidate generator (batch 1: exploration).

Families (all filtered: monic, a0!=0, irreducible over Q, root-count via numpy):
  tri      x^24 + a*x^k + b
  subst    g(x^k) with random deg(24/k) irreducible g
  rand     random dense monic degree 24, small coeffs
  tr       totally-real search (deg 2,3,4,6,8) via random + all-real-roots test
  comp     composita (sum field) of totally-real pairs via resultant
  cheb     T24(x) - c
  cyc      x^24 - a  (a small, 24th-power-free-ish)

Output: one poly per line "a0,...,a24 # family:params r=N"
"""
import random, sys, itertools
import numpy as np
from sympy import Poly, symbols, factor_list, resultant, ZZ

x, y, t = symbols("x y t")

# ---------------- filters ----------------

def count_real_roots(coeffs):
    """coeffs ascending powers; numpy root count (imag tolerance)."""
    r = np.roots(coeffs[::-1])
    return int(np.sum(np.abs(r.imag) < 1e-5))

def is_irreducible(coeffs):
    p = Poly(coeffs, x)
    fl = factor_list(p)
    return len(fl[1]) == 1 and fl[1][0][1] == 1

def to_line(coeffs, tag, r=None):
    if r is None:
        r = count_real_roots(coeffs)
    return ",".join(str(c) for c in coeffs) + f" # {tag} r={r}"

# ---------------- families ----------------

def gen_trinomials(rng, max_a=6, max_b=8, ks=(1,2,3,4,5,6,8,9,10,12,16,18)):
    seen = set()
    for k in ks:
        for a in range(-max_a, max_a+1):
            if a == 0: continue
            for b in range(-max_b, max_b+1):
                if b == 0: continue
                c = [0]*25
                c[0] = b; c[k] = a; c[24] = 1
                key = tuple(c)
                if key in seen: continue
                seen.add(key)
                if not is_irreducible(c): continue
                yield c, f"tri a={a} k={k} b={b}"

def gen_substitution(rng, n_each=60):
    """g(x^k), deg g = d, d*k=24; random irreducible g."""
    for k, d in ((2,12),(3,8),(4,6),(6,4),(8,3),(12,2)):
        made = 0
        tries = 0
        while made < n_each and tries < 2000:
            tries += 1
            g = [0]*(d+1)
            g[d] = 1
            for i in range(d):
                g[i] = rng.randint(-3, 3)
            g[0] = rng.choice([z for z in range(-5,6) if z != 0])
            if not is_irreducible(g): continue
            # f = g(x^k): coeff of x^(k*i) is g[i]
            c = [0]*25
            for i, gi in enumerate(g):
                c[k*i] = gi
            if c[0] == 0: continue
            if not is_irreducible(c): continue
            made += 1
            yield c, f"subst k={k} g={tuple(g)}"

def gen_random_dense(rng, n=120, bound=2):
    made = 0; tries = 0
    while made < n and tries < 6000:
        tries += 1
        c = [rng.randint(-bound, bound) for _ in range(24)] + [1]
        if c[0] == 0: continue
        if not is_irreducible(c): continue
        made += 1
        yield c, "rand_dense"

def gen_totally_real(rng, degrees=(2,3,4,6,8), n_per=40, bound=8):
    """Random search for totally real irreducible polys of small degree."""
    for d in degrees:
        made = 0; tries = 0
        while made < n_per and tries < 20000:
            tries += 1
            c = [rng.randint(-bound, bound) for _ in range(d)] + [1]
            if c[0] == 0: continue
            if not is_irreducible(c): continue
            if count_real_roots(c) != d: continue
            made += 1
            yield c, f"totreal d={d}"

def compositum_sum(fc, gc):
    """Minimal poly (one irreducible factor) of alpha+beta: resultant_y(f(y), g(x-y)).
    Returns list of (coeffs, degree) irreducible factors of degree deg(f)*deg(g)."""
    fy = Poly(fc, y)
    gy = Poly(gc, y)
    gxy = Poly(gc, x - y)
    R = resultant(fy, gxy, y)
    R = Poly(R, x)
    factors = factor_list(R)[1]
    out = []
    for fac, mult in factors:
        cf = fac.all_coeffs()[::-1]  # ascending
        deg = fac.degree()
        if deg == len(fc)-1 + len(gc)-1:
            out.append(cf)
    return out

def gen_composita(tr_pool, rng, n_pairs=120):
    """Composita of totally-real pairs (a,b) with a*b=24, via sum field."""
    # tr_pool: dict deg -> list of coeff lists (totally real, irreducible)
    pairs = [(2,12),(3,8),(4,6),(6,4),(8,3),(12,2)]
    made = 0; tries = 0
    while made < n_pairs and tries < 4000:
        tries += 1
        a, b = rng.choice(pairs)
        pa = tr_pool.get(a); pb = tr_pool.get(b)
        if not pa or not pb: continue
        fa = rng.choice(pa); gb = rng.choice(pb)
        try:
            fs = compositum_sum(fa, gb)
        except Exception:
            continue
        if not fs: continue
        for cf in fs:
            if count_real_roots(cf) != 24: continue
            made += 1
            yield cf, f"comp sum({a}x{b})"
            break

def gen_cheb(n=10):
    # T24(x) coefficients (ascending)
    c = [0]*25
    for j in range(13):
        k = 24 - 2*j
        coef = 2**23 * (24 // (24 - j)) * _binom(24 - j, j)
        c[k] = coef
    c[0] = c[0] - c[0]  # T24 constant term is 0 for even n
    # T24(x) - c0
    out = []
    for c0 in range(2, 2+n):
        cc = c.copy()
        cc[0] = -c0
        if is_irreducible(cc):
            out.append((cc, f"cheb T24-{c0}"))
    return out

def _binom(n, k):
    from math import comb
    return comb(n, k)

def gen_cyclic(rng, n=40):
    made = 0
    for a in range(2, 200):
        if made >= n: break
        if a == 1 or a % 4 == 0: continue
        c = [0]*25; c[0] = -a if rng.random() < 0.5 else a; c[24] = 1
        if not is_irreducible(c): continue
        made += 1
        yield c, f"cyc a={c[0]}"

# ---------------- main ----------------

def main():
    rng = random.Random(20260812)
    out = []
    tags = []

    print("trinomials...", file=sys.stderr)
    for c, tag in gen_trinomials(rng):
        out.append((c, tag))

    print("substitution...", file=sys.stderr)
    for c, tag in gen_substitution(rng):
        out.append((c, tag))

    print("random dense...", file=sys.stderr)
    for c, tag in gen_random_dense(rng):
        out.append((c, tag))

    print("totally-real search...", file=sys.stderr)
    tr_pool = {}
    for c, tag in gen_totally_real(rng):
        d = len(c) - 1
        tr_pool.setdefault(d, []).append(c)
        out.append((c, tag))

    print(f"totally-real pool sizes: { {k: len(v) for k, v in tr_pool.items()} }", file=sys.stderr)

    print("composita...", file=sys.stderr)
    for c, tag in gen_composita(tr_pool, rng, n_pairs=150):
        out.append((c, tag))

    print("chebyshev...", file=sys.stderr)
    for c, tag in gen_cheb(8):
        out.append((c, tag))

    print("cyclic...", file=sys.stderr)
    for c, tag in gen_cyclic(rng):
        out.append((c, tag))

    # dedupe
    seen = set(); final = []
    for c, tag in out:
        key = tuple(c)
        if key in seen: continue
        seen.add(key)
        final.append((c, tag))
    print(f"total candidates: {len(final)}", file=sys.stderr)
    for c, tag in final:
        r = count_real_roots(c)
        print(to_line(c, tag, r))


if __name__ == "__main__":
    main()
