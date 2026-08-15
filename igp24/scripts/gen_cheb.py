#!/usr/bin/env python3
"""Chebyshev / De Moivre substitution: f(T_k(x)) — rotation on the circle.

T_k(cos θ) = cos(kθ): substituting x ↦ T_k(x) is k-fold rotation on the
unit circle (the "circle" structure of the user's conjecture, realized
exactly).  Unlike subst x^k (cyclic wreath), T_k has two ramified points
±1 → the covering group is dihedral-type (D_k), giving different groups.

Configurations (degree 24):
  k=2, base deg 12:  f(T_2(x)) = f(2x²−1)          (12·2)
  k=3, base deg 8:   f(T_3(x)) = f(4x³−3x)         (8·3)
  k=4, base deg 6:   f(T_4(x)) = f(8x⁴−8x²+1)      (6·4)
  k=6, base deg 4:   f(T_6(x)) = f(32x⁶−48x⁴+18x²−1) (4·6)

Bases: cyclotomic-real subfields (deg 12/8/6/4) and quadray fields.
"""
import math, random, sys
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X = pari("x")
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 20261010
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 40
KSUB = int(sys.argv[3]) if len(sys.argv) > 3 else 0  # 0=auto by base deg, else explicit

def is_irred(g):
    try: return bool(g.polisirreducible())
    except Exception: return False

def real_count(g):
    try: return int(pari.polsturm(g))
    except Exception: return -1

def cheb(k):
    """Chebyshev polynomial T_k(x) (exact integer coefficients)."""
    if k == 1:
        return X
    if k == 2:
        return 2 * X**2 - 1
    if k == 3:
        return 4 * X**3 - 3 * X
    if k == 4:
        return 8 * X**4 - 8 * X**2 + 1
    if k == 6:
        return 32 * X**6 - 48 * X**4 + 18 * X**2 - 1
    raise ValueError(f"unsupported k={k}")

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

def monic_ize(c):
    lc = c[-1]
    if lc == 1: return c
    if lc == -1: return [-z for z in c]
    if all(z % lc == 0 for z in c): return [z // lc for z in c]
    d = len(c) - 1
    out = [c[k] * lc ** (d - 1 - k) for k in range(d)]
    out.append(1)
    return out

def compose(f, tk):
    """g(x) = f(T_k(x)) via substitution; root-scaled monic (group-preserving)."""
    g = pari.subst(f, X, tk)
    if g.poldegree() != int(f.poldegree()) * int(tk.poldegree()):
        return None
    c = [int(z) for z in pari.Vecrev(pari.Pol(g, X))]
    c = monic_ize(c)  # root-scaling: field-isomorphic, Galois group preserved
    # root-scaled coeffs can reach 10^110 (T3: 4^(8·23)); allow them — API
    # verification is slower but still works (heavier coeffs seen before)
    if any(abs(z) > 10**130 for z in c):
        return None
    return c

def main():
    rng = random.Random(SEED)
    # base selection: (n, deg) for cyclotomic real subfields
    bases = [(35, 12), (17, 8), (13, 6), (15, 4), (21, 6), (28, 6), (36, 6)]
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 300 * N_TARGET:
        tries += 1
        n, d = rng.choice(bases)
        f0 = get_real_subfield(n, d)
        if f0 is None:
            continue
        k = KSUB if KSUB else {12: 2, 8: 3, 6: 4, 4: 6}[d]
        tk = cheb(k)
        c = compose(f0, tk)
        if c is None:
            continue
        if len(c) != 25:
            continue
        if c[0] == 0:
            continue
        key = tuple(c)
        if key in seen:
            continue
        g = pari.Polrev(c)
        if not is_irred(g):
            continue
        r = int(pari.polsturm(g))
        seen.add(key)
        out.append((c, r, n, k))
        made += 1
    print(f"cheb: made={made} tries={tries}", file=sys.stderr)
    rdist = {}
    kdist = {}
    for c, r, n, k in out:
        rdist[r] = rdist.get(r, 0) + 1
        kdist[f"T{k}(ζ{n}^+)"] = kdist.get(f"T{k}(ζ{n}^+)", 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    print(f"配置分布: {kdist}", file=sys.stderr)
    for c, r, n, k in out:
        print(",".join(str(z) for z in c) + f" # cheb T{k} n={n} r={r}")

if __name__ == "__main__":
    main()
