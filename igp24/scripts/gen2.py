#!/usr/bin/env python3
"""IGP24 candidate generator v2 (cypari2 backend).

Families:
  split    perturbed-split: base = prod(x-r_i)*prod(x^2+p_j) + eps*x^k
           -> irreducible deg d, exactly n_real real roots
  even     f = h(x^2), h totally real w/ positive roots -> deg 2*deg(h) totally real
  special  A3 cubic, S3 cubics, biquadratics, x^4-a
  comp     polcompositum sum-field of deg-a and deg-b bases (a*b=24) -> r = r1*r2
  tri      x^24 + a*x^k + b
"""
import random, sys
from cypari2 import Pari

pari = Pari()
X = pari("x")

def mk(coeffs):
    return pari.Polrev(coeffs)  # ascending powers

def coeffs_of(g):
    # ascending list of integer coeffs (as Python ints)
    return [int(c) for c in pari.Vecrev(pari.Pol(g, X))]

def count_real_roots(g):
    try:
        return int(pari.polsturm(g))
    except Exception:
        return 0

def is_irreducible(g):
    try:
        return bool(g.polisirreducible())
    except Exception:
        return False

def to_line(coeffs, tag, r=None):
    if r is None:
        g = mk(coeffs)
        r = count_real_roots(g)
    return ",".join(str(c) for c in coeffs) + f" # {tag} r={r}"

def mul_linear(poly_list, root):
    out = [0] * (len(poly_list) + 1)
    for i, ci in enumerate(poly_list):
        out[i] -= ci * root
        out[i + 1] += ci
    return out

def mul_quad(poly_list, p_):
    out = [0] * (len(poly_list) + 2)
    for i, ci in enumerate(poly_list):
        out[i] += p_ * ci
        out[i + 2] += ci
    return out

# ---------- perturbed split ----------

def split_base(d, n_real, rng, bound=12, eps_max=6):
    r_i = sorted(rng.sample(range(-bound, bound + 1), n_real))
    if n_real >= 2:
        if min(r_i[i+1] - r_i[i] for i in range(n_real - 1)) < 1:
            return None
    n_pairs = (d - n_real) // 2
    p_j = [rng.randint(1, 9) for _ in range(n_pairs)]
    base = [1]
    for ri in r_i:
        base = mul_linear(base, ri)
    for pj in p_j:
        base = mul_quad(base, pj)
    for _ in range(25):
        k = rng.randint(0, d - 1)
        eps = rng.choice([z for z in range(-eps_max, eps_max + 1) if z != 0])
        c = base.copy()
        c[k] += eps
        if c[0] == 0:
            continue
        g = mk(c)
        if not is_irreducible(g):
            continue
        if count_real_roots(g) == n_real:
            return c
    return None

# ---------- even totally real ----------

def even_totally_real(hdeg, rng, bound=14, eps_max=5):
    d = hdeg
    for _ in range(300):
        r_i = sorted(rng.sample(range(1, bound + 1), d))
        if min(r_i[i+1] - r_i[i] for i in range(d - 1)) < 1:
            continue
        base = [1]
        for ri in r_i:
            base = mul_linear(base, ri)
        for _ in range(25):
            k = rng.randint(0, d - 1)
            eps = rng.choice([z for z in range(-eps_max, eps_max + 1) if z != 0])
            h = base.copy()
            h[k] += eps
            if h[0] == 0:
                continue
            hg = mk(h)
            if count_real_roots(hg) != d:
                continue
            if not is_irreducible(hg):
                continue
            f = [0] * (2 * d + 1)
            for i, ci in enumerate(h):
                f[2 * i] = ci
            fg = mk(f)
            if f[0] == 0 or not is_irreducible(fg):
                continue
            if count_real_roots(fg) != 2 * d:
                continue
            return f
    return None

# ---------- specials ----------

def specials(rng):
    out = []
    for c, tag in (([1, -3, 0, 1], "A3 x^3-3x+1"),
                   ([-1, 3, 0, 1], "A3 -x^3+3x+1")):
        if is_irreducible(mk(c)):
            out.append((c, tag))
    for a in range(1, 9):
        for b in range(1, 9):
            c = [-b, -a, 0, 1]
            if is_irreducible(mk(c)) and count_real_roots(mk(c)) == 3:
                out.append((c, f"S3 x^3-{a}x-{b}"))
    for a in range(-12, -1):
        for b in range(1, 10):
            if a * a <= 4 * b:
                continue
            c = [b, 0, a, 0, 1]
            if is_irreducible(mk(c)):
                out.append((c, f"biq a={a} b={b}"))
    for a in range(2, 8):
        c = [-a, 0, 0, 0, 1]
        if is_irreducible(mk(c)):
            out.append((c, f"x^4-{a}"))
    return out

# ---------- compositum ----------

def compositum_sum(fc, gc):
    f = mk(fc); g = mk(gc)
    try:
        comps = pari.polcompositum(f, g)
    except Exception:
        return []
    deg = len(fc) - 1 + len(gc) - 1
    out = []
    for item in comps:
        R = item[0]
        if R.poldegree() == deg and is_irreducible(R):
            out.append(coeffs_of(R))
    return out

def gen_composita(pools, rng, n_target, want_r=None):
    pairs = [(2,12),(3,8),(4,6),(6,4),(8,3),(12,2)]
    made = 0; tries = 0
    while made < n_target and tries < 4000:
        tries += 1
        a, b = rng.choice(pairs)
        pa = pools.get(a); pb = pools.get(b)
        if not pa or not pb:
            continue
        fa, r1 = rng.choice(pa)
        gb, r2 = rng.choice(pb)
        if want_r is not None and r1 * r2 != want_r:
            continue
        fs = compositum_sum(fa, gb)
        if not fs:
            continue
        for cf in fs:
            r = count_real_roots(mk(cf))
            if want_r is not None and r != want_r:
                continue
            made += 1
            yield cf, f"comp {a}x{b} ({r1}*{r2})"
            break

# ---------- main ----------

def main():
    rng = random.Random(20260813)
    out = []
    pools = {}

    print("building pools...", file=sys.stderr)
    for d, n_real in ((2,2),(3,3),(3,1),(4,4),(4,2),(4,0),(6,6),(6,4),(6,2),
                      (8,8),(8,6),(8,4),(8,2),(12,12),(12,8),(12,4),(12,0)):
        got = 0; tries = 0
        while got < 25 and tries < 6000:
            tries += 1
            c = split_base(d, n_real, rng)
            if c is None:
                continue
            pools.setdefault(d, []).append((c, n_real))
            got += 1
        print(f"  split d={d} n_real={n_real}: {got}", file=sys.stderr)

    for hd in (3, 4, 6):
        got = 0; tries = 0
        while got < 20 and tries < 600:
            tries += 1
            f = even_totally_real(hd, rng)
            if f is None:
                continue
            pools.setdefault(2*hd, []).append((f, 2*hd))
            got += 1
        print(f"  even d={2*hd}: {got}", file=sys.stderr)

    for c, tag in specials(rng):
        d = len(c) - 1
        pools.setdefault(d, []).append((c, count_real_roots(mk(c))))
        out.append((c, tag))

    for k, v in pools.items():
        print(f"pool deg {k}: {len(v)}", file=sys.stderr)

    for want_r, n in ((24, 260), (16, 160), (20, 120), (12, 80), (8, 60),
                      (4, 40), (0, 40), (18, 30), (6, 30)):
        print(f"composita r={want_r} x{n}...", file=sys.stderr)
        for c, tag in gen_composita(pools, rng, n, want_r=want_r):
            out.append((c, tag))

    print("trinomials...", file=sys.stderr)
    for k in (1,2,3,4,5,6,8,9,10,12,16,18):
        for a in (-5,-4,-3,-2,2,3,4,5):
            for b in (-6,-5,-4,-3,-2,2,3,4,5,6):
                c = [0]*25; c[0]=b; c[k]=a; c[24]=1
                if is_irreducible(mk(c)):
                    out.append((c, f"tri a={a} k={k} b={b}"))

    seen = set(); final = []
    for c, tag in out:
        key = tuple(c)
        if key in seen:
            continue
        seen.add(key)
        final.append((c, tag))
    print(f"total: {len(final)}", file=sys.stderr)
    for c, tag in final:
        print(to_line(c, tag))

if __name__ == "__main__":
    main()
