#!/usr/bin/env python3
"""升维塔生成器 B: 全实 abelian degree-12 基域 + √ 层 (degree 24).

数学 (2026-08-13 分析):
- 开放对 65% 的 2-部分 ∈ [2^10, 2^14] (3,177 对), 73% 只要 3-部分 3^1
- helix 家族用复数分圆基域 (Q(ζ3), degree 2) → 2-部分锁死在 2^5-2^7, r 锁死 0
- 升维 = 基域换全实 abelian degree-12 分圆子域 (Q(ζ35)^+, Gal order 12):
  * √ 层元素的共轭数 = 基域嵌入数 = 12 → 闭包 2-部分 ≈ 2^12·2^2 = 2^14 (靶区)
  * 基域全实 (12 实嵌入) → r = 2·#(σ(b) > 0), b 的符号可控
  * b = γ²+c (全正共轭) → r=24; 随机 b → r≈12 (二项); γ²-c → r 小
- 用户直觉: 圆 (分圆) 的正确用法是实子域 = 圆的上半弧; 旋转群 C12 比 C2 大
- 3-部分: 基域 abelian order 12 = 2^2·3^1 → 3^1 (够 73% 开放对)

用法: python3 -u gen_liftB.py <seed> <target> <mode> [n] [bound]
  mode: full (γ²+c → r=24) | mix (随机 b) | neg (γ²-c → r 小)
"""
import math, random, sys, time
sys.set_int_max_str_digits(10**7)
from cypari2 import Pari

pari = Pari()
X = pari("x")
Y = pari("y")
Z = pari("z")

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
N_TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 30
MODE = sys.argv[3] if len(sys.argv) > 3 else "full"
N = int(sys.argv[4]) if len(sys.argv) > 4 else 35
BOUND = int(sys.argv[5]) if len(sys.argv) > 5 else 2


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


def get_real_subfield(n, deg):
    """polsubcyclo(n, deg) 中全实的那个子域 (sturm=deg)."""
    res = pari.polsubcyclo(n, deg)
    polys = res if res.type() == "t_VEC" else [res]
    for p in polys:
        if p.poldegree() == deg and int(pari.polsturm(p)) == deg:
            return p
    return None


def monic_ize(c):
    """整系数向量 → 首一多项式系数 (a0..a_{d-1}, a_d=1)."""
    lc = c[-1]
    if lc == 1:
        return c
    if lc == -1:
        return [-z for z in c]
    if all(z % lc == 0 for z in c):
        return [z // lc for z in c]
    d = len(c) - 1
    out = [c[k] * lc ** (d - 1 - k) for k in range(d)]
    out.append(1)
    return out


def sqrt_layer(f, bcoeffs):
    """Res_y(f(y), z² - b(y)) → degree 2·deg(f) 首一多项式系数 (a0..).

    注意: bcoeffs 是升幂系数 (a0 + a1·y + ...), 而 PARI Pol() 收降幂
    → 必须反转: Pol(bcoeffs[::-1], Y).
    """
    d = f.poldegree()
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
    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
    c = prim_part(c)
    P = pari.Polrev(c)
    if not is_irred(P):
        return None
    cc = [int(z) for z in pari.Vecrev(pari.Pol(P, X))]
    if any(abs(z) > 10**40 for z in cc):
        return None
    return monic_ize(cc)


def b_elem(f, rng, cmin, cmax, deg=2):
    """低次小系数元素 + 常数偏移: b = Σ bc[j]·x^j, bc[j]∈{-1,0,1}, 偏移 ∈ [cmin,cmax].

    σ(b) = Σ bc[j]·α^j + 偏移, |α| ≤ 2 → 非偏移部分 ≤ 2^(deg+1)-2.
    偏移控制符号 → 控制 r. 系数小 → resultant 系数 ~10^14-10^18 (可提交).
    """
    d = int(f.poldegree())
    for _ in range(60):
        bc = [0] * d
        k = rng.randint(0, min(deg, d - 1))
        for j in range(1, k + 1):
            bc[j] = rng.randint(-1, 1)
        if all(z == 0 for z in bc):
            bc[0] = 1
        bc[0] += rng.randint(cmin, cmax)
        return bc
    return None


def b_full_positive(f, rng):
    """偏移 8..12 → σ(b) ≥ 1 全正 → r = 24."""
    return b_elem(f, rng, 8, 12)


def b_high(f, rng):
    """偏移 0..3 → 8-11 个正共轭 → r ∈ {14,16,18,20} (r=16/20 缺口)."""
    return b_elem(f, rng, 0, 3)


def b_random(f, rng):
    """偏移 -4..8 → r 多样."""
    return b_elem(f, rng, -4, 8)


def b_negative(f, rng):
    """偏移 -12..-8 → σ(b) < 0 全负 → r = 0."""
    return b_elem(f, rng, -12, -8)


def main():
    rng = random.Random(SEED)
    t0 = time.time()
    f0 = get_real_subfield(N, 12)
    if f0 is None:
        print(f"FAIL: no real degree-12 subfield of Q(ζ{N})", file=sys.stderr)
        return
    print(f"f0 = real subfield Q(ζ{N})^+, deg={f0.poldegree()} sturm={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 200 * N_TARGET:
        tries += 1
        if MODE == "full":
            bc = b_full_positive(f0, rng)
        elif MODE == "hi":
            bc = b_high(f0, rng)
        elif MODE == "neg":
            bc = b_negative(f0, rng)
        else:
            bc = b_random(f0, rng)
        if bc is None:
            continue
        c = sqrt_layer(f0, bc)
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
        out.append((c, r))
        made += 1
    dt = time.time() - t0
    print(f"liftB_{MODE}_n{N}: made={made} tries={tries} {dt:.0f}s", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # liftB_{MODE}_n{N} r={r}")


if __name__ == "__main__":
    main()
