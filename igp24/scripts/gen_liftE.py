#!/usr/bin/env python3
"""升维塔 E: Q(ζ13)^+(√b)(√c) — 两层 √ 共轭爆炸 (degree 6·2·2 = 24).

数学 (2026-08-13):
- liftB (单层 √, 基域 order-12 abelian) 产出 24T14744/14745 全被覆盖
- 塔 E 基域 Q(ζ13)^+ = Q(2cos(2π/13)), degree 6, Gal = C6 (2^1·3^1)
- 第一层 √b: b ∈ K0 共轭 6 → 闭包 2^6
- 第二层 √c: c ∈ K1 (degree 12) 共轭 12 → 闭包 2^12
- 总 2-部分 ≈ 2^6·2^12·2^1 = 2^19 → 打 2^16·3^1 (359 对) / 2^17·3^1 (174) / 2^18·3^1 (121) 桶
- r 控制: 每层 √ 的符号 → 两层全正偏移 → r=24 (全实缺口 1,221 对)
- c 用 γ²+d (d∈{1,2,3,5}) 保证全正; 随机 d 或负偏移给 r 多样性

用法: python3 -u gen_liftE.py <seed> <target> <mode> [bound]
  mode: full (两层全正 → r=24) | mix (随机) | neg (r=0)
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
BOUND = int(sys.argv[4]) if len(sys.argv) > 4 else 2
BASE = sys.argv[5] if len(sys.argv) > 5 else "cyclo13"  # cyclo13 | s3x2


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
    res = pari.polsubcyclo(n, deg)
    polys = res if res.type() == "t_VEC" else [res]
    for p in polys:
        if p.poldegree() == deg and int(pari.polsturm(p)) == deg:
            # polredabs 压系数 → 提高 resultant 产率
            try:
                pr = pari.polredabs(p)
                if pr.poldegree() == deg and int(pari.polsturm(pr)) == deg:
                    return pr
            except Exception:
                pass
            return p
    return None


def get_s3xc2_base(d):
    """Q(√d)(α), α = x³-4x+1 的根 (全实 S3, disc=229 非平方), degree 6, Gal = S3×C2 (非 abelian)."""
    f_alpha = pari.Polrev([1, -4, 0, 1])  # 升幂: 1 - 4x + x³, 全实 S3
    if int(pari.polsturm(f_alpha)) != 3:
        return None
    comp = pari.polcompositum(f_alpha, X**2 - d)
    for p in comp:
        if p.poldegree() == 6:
            try:
                pr = pari.polredabs(p)
                if pr.poldegree() == 6:
                    return pr
            except Exception:
                return p
    return None


def get_c9xc2_base(d):
    """Q(ζ9)^+(√d): Q(ζ9)^+ = x³-3x+1 (全实 C3), 二次扩张 degree 6, 全实, Gal=C6 (域不同 → 群不同)."""
    f_zeta9 = pari.Polrev([1, -3, 0, 1])  # x³ - 3x + 1, 全实 C3
    if int(pari.polsturm(f_zeta9)) != 3:
        return None
    comp = pari.polcompositum(f_zeta9, X**2 - d)
    for p in comp:
        if p.poldegree() == 6:
            try:
                pr = pari.polredabs(p)
                if pr.poldegree() == 6:
                    return pr
            except Exception:
                return p
    return None


def monic_ize(c):
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


def sqrt_layer(f, bcoeffs, do_polred=False):
    """Res_y(f(y), z² - b(y)) → degree 2·deg(f) 系数 (a0..). bcoeffs 升幂."""
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
    c = [int(z) for z in pari.Vecrev(pari.Pol(R, X))]
    c = prim_part(c)
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


def b_lowdeg(f, rng, cmin, cmax, deg=2):
    """低次小系数 + 偏移. 升幂系数."""
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


def b_square_shift(f, rng, dvals=(1, 2, 3, 5), sign=1):
    """b = γ² + d (或 -γ²-d): 全正/全负共轭. γ 必须低次 (deg≤2) 否则共轭值爆炸."""
    dd = int(f.poldegree())
    for _ in range(60):
        gc = [rng.randint(-1, 1) for _ in range(min(3, dd))] + [0] * (dd - min(3, dd))
        if all(z == 0 for z in gc):
            continue
        gam = pari.Polrev(gc, X)
        b = (gam * gam + rng.choice(dvals)) % f
        bc = [int(z) for z in pari.Vecrev(pari.Pol(b, X))]
        while len(bc) < dd:
            bc.append(0)
        if sign < 0:
            bc = [-z for z in bc]
        return bc
    return None


def main():
    rng = random.Random(SEED)
    t0 = time.time()
    if BASE == "cyclo13":
        f0 = get_real_subfield(13, 6)
        tag_base = "Q(ζ13)^+"
    elif BASE == "c9xc2":
        f0 = get_c9xc2_base(rng.choice([2, 3, 5, 6, 7]))
        tag_base = "Q(ζ9)^+(√d)"
    else:
        f0 = get_s3xc2_base(rng.choice([2, 3, 5, 6, 7]))
        tag_base = "S3×C2"
    if f0 is None:
        print(f"FAIL: no base for {BASE}", file=sys.stderr)
        return
    print(f"f0 = {tag_base} deg={f0.poldegree()} sturm={int(pari.polsturm(f0))}", file=sys.stderr)
    out = []
    seen = set()
    made = 0
    tries = 0
    while made < N_TARGET and tries < 400 * N_TARGET:
        tries += 1
        # L1 = L0(√b)
        if MODE == "full":
            # b_lowdeg 偏移 8-12: σ(b) = 偏移 + Σbc[j]σ(x)^j, |σ(x)|≤2, deg≤2 → σ(b)∈(2,14) 全正
            # 系数小 (10^9) → nfdisc 快; γ²+d 系数 10^26 太慢, 弃用
            bc = b_lowdeg(f0, rng, 8, 12)
        elif MODE == "neg":
            bc = b_lowdeg(f0, rng, -12, -8)  # 全负 → r=0
        else:
            bc = b_lowdeg(f0, rng, -4, 8)  # 随机
        if bc is None:
            continue
        c1 = sqrt_layer(f0, bc, do_polred=True)
        if c1 is None:
            continue
        f1 = pari.Polrev(c1)
        # L2 = L1(√c)
        if MODE == "full":
            cc = b_lowdeg(f1, rng, 8, 12)  # K1 根尺度可能>2, 偏移 8-12 需验证全正
        elif MODE == "neg":
            cc = b_lowdeg(f1, rng, -12, -8)
        else:
            cc = b_lowdeg(f1, rng, -4, 8)
        if cc is None:
            continue
        c2 = sqrt_layer(f1, cc)
        if c2 is None:
            continue
        if len(c2) != 25:
            continue
        if c2[0] == 0:
            continue
        key = tuple(c2)
        if key in seen:
            continue
        g = pari.Polrev(c2)
        if not is_irred(g):
            continue
        r = int(pari.polsturm(g))
        seen.add(key)
        out.append((c2, r))
        made += 1
    dt = time.time() - t0
    print(f"liftE_{MODE}: made={made} tries={tries} {dt:.0f}s", file=sys.stderr)
    rdist = {}
    for c, r in out:
        rdist[r] = rdist.get(r, 0) + 1
    print(f"r分布: {dict(sorted(rdist.items()))}", file=sys.stderr)
    for c, r in out:
        print(",".join(str(z) for z in c) + f" # liftE_{MODE} r={r}")


if __name__ == "__main__":
    main()
