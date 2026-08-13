#!/usr/bin/env python3
"""升维塔可行性验证 (2026-08-13)

验证三个"升维"塔结构是否产出 degree-24 不可约多项式, 以及 r 分布:
  A: Q(ζ13)(√b)          — degree 12·2=24, 基域 Gal=C12 (2^2·3), √ 层共轭爆炸
  B: Q(ζ35)^+(√b)        — degree 12·2=24, 全实 abelian 基域, r 缺口靶
  C: Q(ζ5)(∛a)(√b)       — degree 4·3·2=24, 旋转层 C4 + Kummer + √ 爆炸
"""
import math, random, sys, time
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X = pari("x")
Y = pari("y")
Z = pari("z")

def is_irred(g):
    try:
        return bool(g.polisirreducible())
    except Exception:
        return False

def prim_part(c):
    gg = 0
    for z in c:
        gg = math.gcd(gg, abs(z))
    if gg > 1:
        return [z // gg for z in c]
    return c

def poly_layer(f, bcoeffs, kth, do_polred=False):
    """Res_y(f(y), z^k - b(y)) with z→x; monic degree k*deg(f) poly."""
    d = f.poldegree()
    fy = pari.subst(f, X, Y)
    b = pari.Pol(bcoeffs, Y)
    coeffs = [Z**kth - bcoeffs[0]] + [-bcoeffs[i] for i in range(1, d)]
    g = pari.Pol(coeffs, Y)
    t0 = time.time()
    try:
        R = pari.polresultant(fy, g, Y)
    except Exception as e:
        print("  resultant fail:", e, file=sys.stderr)
        return None
    dt = time.time() - t0
    R = pari.subst(R, Z, X)
    if R.poldegree() != kth * d:
        print(f"  deg mismatch {R.poldegree()} != {kth*d}", file=sys.stderr)
        return None
    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
    c = prim_part(c)
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    if do_polred and kth * d <= 12:
        try:
            Pr = pari.polredabs(P)
            if Pr.poldegree() == kth * d:
                P = Pr
        except Exception:
            pass
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return cc

def rand_elem(f, rng, bound=4):
    d = f.poldegree()
    for _ in range(40):
        bc = [rng.randint(-bound, bound) for _ in range(d)]
        if all(z == 0 for z in bc):
            continue
        return bc
    return None

def tower_A(seed, n_target=8):
    """Q(ζ13)(√b): f0 = Φ13(x) degree 12, then sqrt layer."""
    rng = random.Random(seed)
    f0 = pari.polcyclo(13)
    out = []
    for _ in range(30):
        bc = rand_elem(f0, rng, bound=2)
        if bc is None:
            continue
        c = poly_layer(f0, bc, 2)
        if c is None:
            continue
        if len(c) != 25:
            continue
        g = pari.Polrev(c)
        r = int(pari.polsturm(g))
        out.append((c, r))
        if len(out) >= n_target:
            break
    print(f"[A Q(ζ13)(√b)] made={len(out)}")
    for c, r in out[:n_target]:
        print(f"A r={r} len={len(c)} first3={c[:3]}")

def tower_B(seed, n_target=8):
    """Q(ζ35)^+(√b): real subfield of Q(ζ35), degree 12, then sqrt."""
    rng = random.Random(seed)
    # polsubcyclo(35, 12)?  Q(ζ35) 的 degree-12 实子域 = 最大实子域 (ζ35+ζ35^-1)
    # 用 polsubcyclo(n, d): n=35, 找 degree-12 子域
    res = pari.polsubcyclo(35, 12)
    polys = res if res.type() == "t_VEC" else [res]
    f0 = None
    for p in polys:
        if p.poldegree() == 12:
            f0 = p
            break
    if f0 is None:
        print("[B] no degree-12 subfield of Q(ζ35) found")
        return
    print(f"[B] f0 deg={f0.poldegree()} sturm={int(pari.polsturm(f0))} (全实={int(pari.polsturm(f0))==12})")
    out = []
    for _ in range(30):
        bc = rand_elem(f0, rng, bound=2)
        if bc is None:
            continue
        c = poly_layer(f0, bc, 2)
        if c is None:
            continue
        if len(c) != 25:
            continue
        g = pari.Polrev(c)
        r = int(pari.polsturm(g))
        out.append((c, r))
        if len(out) >= n_target:
            break
    print(f"[B Q(ζ35)^+(√b)] made={len(out)} r分布={sorted(set(r for _, r in out))}")
    for c, r in out[:n_target]:
        print(f"B r={r} len={len(c)}")

def tower_C(seed, n_target=8):
    """Q(ζ5)(∛a)(√b): f0 = Φ5 degree 4, cbrt layer (ζ3 ∉ Q(ζ5), 闭包要 ζ3), sqrt layer."""
    rng = random.Random(seed)
    f0 = pari.polcyclo(5)
    out = []
    for _ in range(20):
        a = rng.randint(2, 12)
        c1 = poly_layer(f0, [0] * 4, 3) if False else None
        # 直接用 a ∈ Q: x³ - a, resultant with f0
        fy = pari.subst(f0, X, Y)
        g = Z**3 - a
        try:
            R = pari.polresultant(fy, pari.subst(g, Z, Y), Y)
        except Exception:
            continue
        R = pari.subst(R, Z, X)
        if R.poldegree() != 12:
            continue
        cc = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
        cc = prim_part(cc)
        f1 = pari.Polrev(cc)
        if not is_irred(f1):
            continue
        # L2 = L1(√b), b random in L1
        bc = rand_elem(f1, rng, bound=2)
        if bc is None:
            continue
        c = poly_layer(f1, bc, 2)
        if c is None:
            continue
        if len(c) != 25:
            continue
        g = pari.Polrev(c)
        r = int(pari.polsturm(g))
        out.append((c, r))
        if len(out) >= n_target:
            break
    print(f"[C Q(ζ5)(∛a)(√b)] made={len(out)} r分布={sorted(set(r for _, r in out))}")
    for c, r in out[:n_target]:
        print(f"C r={r} len={len(c)}")

if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    which = sys.argv[2] if len(sys.argv) > 2 else "ABC"
    t0 = time.time()
    if "A" in which:
        tower_A(seed)
        print(f"--- A done {time.time()-t0:.0f}s")
    if "B" in which:
        tower_B(seed)
        print(f"--- B done {time.time()-t0:.0f}s")
    if "C" in which:
        tower_C(seed)
        print(f"--- C done {time.time()-t0:.0f}s")
