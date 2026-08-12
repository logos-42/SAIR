#!/usr/bin/env python3
"""IGP24 candidate generator v3: cyclotomic subfields (abelian groups).

polsubcyclo(n, d) gives a degree-d subfield of Q(zeta_n). For primes p with
24 | (p-1)/2 and -1 in the relevant subgroup, polsubcyclo(p, 24) is a totally
real cyclic C24 field (r=24). Composite n give non-cyclic abelian groups.
Composita of cyclotomic subfields give more abelian group variety.
"""
import random, sys
from cypari2 import Pari

pari = Pari()
X = pari("x")

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True

def poly_ok(f):
    try:
        return f.poldegree() == 24 and f.polisirreducible()
    except Exception:
        return False

def subcyclo24(n):
    """All degree-24 subfields of Q(zeta_n); return (poly, real_roots) list."""
    try:
        res = pari.polsubcyclo(n, 24)
    except Exception:
        return []
    if res.type() == "t_VEC":
        items = list(res)
    else:
        items = [res]
    out = []
    for f in items:
        if not poly_ok(f):
            continue
        r = int(pari.polsturm(f))
        out.append((f, r))
    return out

def main():
    rng = random.Random(20260814)
    out = []
    seen = set()

    # 1. primes p = 1 mod 48: cyclic C24 totally real (r=24)
    primes24 = []
    p = 97
    while p < 4000:
        if is_prime(p) and (p - 1) % 48 == 0:
            primes24.append(p)
        p += 48
    print(f"primes p≡1 mod 48 < 4000: {len(primes24)}", file=sys.stderr)
    for p in primes24:
        for f, r in subcyclo24(p):
            c = [int(z) for z in pari.Vecrev(pari.Pol(f, X))]
            key = tuple(c)
            if key in seen: continue
            seen.add(key)
            out.append((c, f"cyclo C24 p={p}"))

    # 2. primes p = 1 mod 24 (but not 48): C24 with r=0 or r=2 etc.
    p = 73
    n = 0
    while p < 3000 and n < 40:
        if is_prime(p) and (p - 1) % 24 == 0 and (p - 1) % 48 != 0:
            for f, r in subcyclo24(p):
                c = [int(z) for z in pari.Vecrev(pari.Pol(f, X))]
                key = tuple(c)
                if key in seen: continue
                seen.add(key)
                out.append((c, f"cyclo C24b p={p}"))
                n += 1
        p += 24

    # 3. composite n: non-cyclic abelian groups
    cand = []
    for a in range(3, 40):
        for b in range(3, 40):
            n = a * b
            if n <= 24 or n > 1500: continue
            if is_prime(a) and is_prime(b):
                cand.append(n)
    cand = sorted(set(cand))[:80]
    for n in cand:
        try:
            phi = pari.eulerphi(n)
            if phi % 24 != 0: continue
        except Exception:
            continue
        for f, r in subcyclo24(n):
            c = [int(z) for z in pari.Vecrev(pari.Pol(f, X))]
            key = tuple(c)
            if key in seen: continue
            seen.add(key)
            out.append((c, f"cyclo ab n={n}"))

    # 4. cyclotomic degree-12 real subfields (for C12 x C2 composita)
    pool12 = []
    p = 73
    while p < 2000 and len(pool12) < 30:
        if is_prime(p) and (p - 1) % 24 == 0:
            try:
                res = pari.polsubcyclo(p, 12)
            except Exception:
                p += 24; continue
            items = list(res) if res.type() == "t_VEC" else [res]
            for f in items:
                try:
                    if f.poldegree() == 12 and f.polisirreducible() and int(pari.polsturm(f)) == 12:
                        pool12.append([int(z) for z in pari.Vecrev(pari.Pol(f, X))])
                except Exception:
                    continue
        p += 24
    print(f"C12 real pool: {len(pool12)}", file=sys.stderr)
    # quadratics x^2 - a, a squarefree positive
    quads = []
    for a in range(2, 60):
        if a % 4 == 0: continue
        fa = pari.Polrev([-a, 0, 1])
        if fa.polisirreducible():
            quads.append([-a, 0, 1])
    # composita C12 x C2 (totally real)
    made = 0; tries = 0
    while made < 60 and tries < 2000:
        tries += 1
        g12 = rng.choice(pool12)
        q2 = rng.choice(quads)
        try:
            comps = pari.polcompositum(pari.Polrev(g12), pari.Polrev(q2))
        except Exception:
            continue
        for R in comps:
            try:
                if R.poldegree() == 24 and R.polisirreducible() and int(pari.polsturm(R)) == 24:
                    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
                    key = tuple(c)
                    if key in seen: continue
                    seen.add(key)
                    out.append((c, "cyclo C12xC2"))
                    made += 1
                    break
            except Exception:
                continue

    print(f"total gen3: {len(out)}", file=sys.stderr)
    for c, tag in out:
        r = int(pari.polsturm(pari.Polrev(c)))
        print(",".join(str(z) for z in c) + f" # {tag} r={r}")

if __name__ == "__main__":
    main()
