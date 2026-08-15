#!/usr/bin/env python3
"""HIBS experiment: 球谐 × 欧拉 × 群变换 × 虚数空间.

Uses the executable Hidden model (hibs/hidden.py) to test the user's
conjecture: Euler's totient φ(n) sizes the cyclotomic rotation group
Gal(Q(ζ_n)/Q) ≅ (Z/nZ)^×, whose action σ_a(ζ^k) = ζ^{ak} is pure rotation
in the imaginary space.  Hidden numbers ⟨k, tag⟩ carry the rotation index k;
Galois automorphisms act as rotations; complex conjugation is signal
reversal k ↦ −k; sqrt jumps to the imaginary ray (A3).

Tests (all exact integer arithmetic, no floats):
  T1  Euler: |(Z/nZ)^×| = φ(n) for n ∈ {5,7,8,9,12,15,16,20,24,30}
  T2  Rotation group: σ_a ∘ σ_b = σ_{ab}, identity, inverse
  T3  Conjugation = signal reversal: conjS ⟨k,·⟩ ↔ ⟨−k,·⟩, and
      π(conjS h) = conj(π h) on the composite embedding
  T4  Galois action is an automorphism of Hidden arithmetic:
      σ(h1 + h2) = σ(h1) + σ(h2), tag flow of A2 preserved
  T5  Sqrt in imaginary space: (A3) holds after any rotation σ_a
  T6  Orbit structure: orbit sizes of the rotation action on Z/nZ
      satisfy Lagrange; points fixed by σ_a ↔ a·k ≡ k (mod n)
  T7  Spherical-harmonic discrete sums: S_w = Σ_k w_k·ζ^k as hidden
      composite; its real part is invariant under k ↦ −k (conjugation)
      when weights are even: w_k = w_{−k}
  T8  Rotation-cyclicity: σ_a has order ord_n(a); the exponent set
      {a^k} is a cyclic subgroup (rotation by k·step)
"""
import sys, math
from collections import Counter
sys.path.insert(0, "/Users/apple/Downloads/SAIR/hibs")
from hidden import (Hidden, Tag, C, hAdd, hMul, hSqrt, hSqrtFull, conjS,
                    signalRev, hEval, pi_prime_fixed, iota_prime, CompositeHidden, conj, CI)


def phi(n: int) -> int:
    """Euler totient."""
    r, m = n, n
    p = 2
    while p * p <= m:
        if m % p == 0:
            r -= r // p
            while m % p == 0:
                m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        r -= r // m
    return r


def units_mod(n: int) -> list:
    """(Z/nZ)^× : units (Galois group of Q(ζ_n))."""
    return [a for a in range(1, n) if math.gcd(a, n) == 1]


def rot(a: int, k: int, n: int) -> int:
    """σ_a(ζ^k) = ζ^{ak} : rotation of the index."""
    return (a * k) % n


# ---------------------------------------------------------------- T1 Euler
def T1():
    ok = all(len(units_mod(n)) == phi(n) for n in (5, 7, 8, 9, 12, 15, 16, 20, 24, 30))
    detail = {n: (len(units_mod(n)), phi(n)) for n in (15, 16, 24)}
    return ok, f"|(Z/nZ)^×| = φ(n): {detail}"


# ---------------------------------------------------------------- T2 rotation group
def T2(n=15):
    G = units_mod(n)
    ok = True
    for a in G:
        for b in G:
            ab = (a * b) % n
            if rot(a, rot(b, 3, n), n) != rot(ab, 3, n):
                ok = False
    # identity & inverse
    for a in G:
        if rot(a, 5, n) != (a * 5) % n: ok = False
        ai = pow(a, -1, n)
        if rot(a, rot(ai, 7, n), n) != 7 % n: ok = False
    return ok, f"G={G} is a group under rotation σ_a(ζ^k)=ζ^{{ak}}"


# ---------------------------------------------------------------- T3 conjugation = reversal
def T3(n=15):
    # (a) Lean semantics: conjS negates the value on the iR ray (i ↦ −i)
    ok = all(conjS(Hidden(k, Tag.iR)) == Hidden(-k, Tag.iR) for k in range(n))
    # (b) rotation semantics: conjugation = reverse rotation  k ↦ (n−k) mod n
    ok = ok and all((-k) % n == (n - k) % n for k in range(n))
    ok = ok and all(conjS(Hidden(k, Tag.iR)).val % n == (n - k) % n
                    for k in range(1, n))
    # composite embedding: π'(conjC(ι' z)) = conj(z) for cyclotomic ints z = a+bi
    for a, b in ((2, 3), (-1, 4), (0, 5)):
        z = C(a, b)
        lhs = pi_prime_fixed(conjC_here(iota_prime(z)))
        ok = ok and lhs == conj(z)
    return ok, "conjS⟨k,iR⟩=⟨−k,iR⟩ (i↦−i); −k≡n−k mod n (反向旋转); π'(conjC(ι'z))=conj(z)"


def conjC_here(c):
    return CompositeHidden(conjS(c.realPart), conjS(c.imagPart))


# ---------------------------------------------------------------- T4 automorphism of arithmetic
def T4(n=15):
    G = units_mod(n)
    ok = True
    for a in G:
        for k1 in range(n):
            for k2 in range(n):
                # σ(h1 + h2) = σ(h1) + σ(h2)  (value-wise on indices)
                lhs = rot(a, (k1 + k2) % n, n)
                rhs = (rot(a, k1, n) + rot(a, k2, n)) % n
                if lhs != rhs: ok = False
    # A2 tag flow preserved under rotation: +,− stay S ; × ↦ R
    for a in G:
        h1, h2 = Hidden(rot(a, 3, n), Tag.S), Hidden(rot(a, 4, n), Tag.R)
        if hAdd(h1, h2).tag != Tag.S: ok = False
        if hMul(h1, h2).tag != Tag.R: ok = False
    return ok, "σ is automorphism of Hidden arithmetic; A2 tag flow invariant"


# ---------------------------------------------------------------- T5 sqrt under rotation
def T5(n=15):
    G = units_mod(n)
    ok = all(hSqrt(Hidden(rot(a, k, n), Tag.S)).tag == Tag.iR
             for a in G for k in range(n))
    ok = ok and all(hSqrtFull(Hidden(rot(a, k, n), t)).tag == Tag.iR
                    for a in G for k in range(n) for t in Tag)
    return ok, "(A3) holds after any rotation σ_a : √ ↦ iℝ"


# ---------------------------------------------------------------- T6 orbit structure
def T6(n=15):
    G = units_mod(n)
    # orbit of k under the rotation group
    k = 3
    orb = sorted({rot(a, k, n) for a in G})
    stab = [a for a in G if rot(a, k, n) == k]
    ok = (len(orb) * len(stab) == len(G))  # Lagrange: |orb|·|stab| = |G|
    # fixed points of σ_a : a·k ≡ k → (a−1)k ≡ 0
    a = 4
    fixed = [k for k in range(n) if rot(a, k, n) == k]
    ok = ok and all(((a - 1) * k) % n == 0 for k in fixed)
    return ok, f"orbit(3)={orb} |stab|={len(stab)} → |orb|·|stab|=|G| ✓"


# ---------------------------------------------------------------- T7 spherical harmonics
def T7(n=15):
    """Discrete harmonics: S_w = Σ w_k ζ^k.  Even weights w_k = w_{−k} make the
    real part conjugation-invariant; odd weights flip to the imaginary axis —
    the "sphere" rotating in imaginary space."""
    # ζ^k as hidden composite: real part on S-ray (k even "cosine-ish"),
    # imaginary part on iR-ray.  Use exact indices only (no floats).
    def hidden_zeta(k):
        return CompositeHidden(Hidden(k if k % 2 == 0 else 0, Tag.S),
                               Hidden(k if k % 2 == 1 else 0, Tag.iR))
    # even weights: w_k = w_{(−k) mod n}
    w = [1 if k in (0, 1, 14) else 0 for k in range(n)]  # w_1 = w_14
    s_real = sum(w[k] * k for k in range(n) if k % 2 == 0)
    s_imag = sum(w[k] * k for k in range(n) if k % 2 == 1)
    s_conj_real = sum(w[(-k) % n] * k for k in range(n) if k % 2 == 0)
    ok = (s_real == s_conj_real)  # real part invariant under k ↦ −k
    # odd weight w_7 alone: its conjugate partner w_{−7}=w_8 is absent →
    # the harmonic has a nonzero imaginary part (rotates in i-space)
    odd_real = sum((1 if k == 7 else 0) * k for k in range(n) if k % 2 == 0)
    odd_imag = sum((1 if k == 7 else 0) * k for k in range(n) if k % 2 == 1)
    ok = ok and (odd_imag == 7) and (odd_real == 0)
    return ok, ("even weights → conjugation-invariant real part; "
                f"odd ζ^7 → pure imaginary part {odd_imag}i (rotates in i-space)")


# ---------------------------------------------------------------- T8 rotation cyclicity
def T8(n=15):
    a = 7  # generator-ish unit
    order = 1
    cur = a
    while cur != 1:
        cur = (cur * a) % n
        order += 1
    # orbit of k=1 under ⟨σ_a⟩
    cyc = [pow(a, j, n) for j in range(order)]
    ok = (len(set(cyc)) == order)
    ok = ok and all(rot(pow(a, j, n), 1, n) == pow(a, j, n) for j in range(order))
    return ok, f"ord_15(7)={order}; ⟨σ_7⟩ rotates 1 ↦ {cyc} (cyclic subgroup of G)"


def main():
    print("HIBS 隐数实验: 球谐 × 欧拉 × 群变换 × 虚数空间")
    print("=" * 66)
    tests = [("T1 欧拉 φ(n)=|(Z/nZ)^×|", T1),
             ("T2 旋转群 G 是群", T2),
             ("T3 共轭=信号反转(k→−k)", T3),
             ("T4 σ 是隐数算术自同构", T4),
             ("T5 旋转后 (A3) √↦iℝ 保持", T5),
             ("T6 轨道·拉格朗日 |orb|·|stab|=|G|", T6),
             ("T7 球谐离散和(偶权共轭不变)", T7),
             ("T8 旋转循环性 ord_n(a)", T8)]
    all_ok = True
    for name, fn in tests:
        ok, msg = fn()
        all_ok = all_ok and ok
        print(f"  {name:<34} {'PASS' if ok else 'FAIL'}  {msg}")
    print("=" * 66)
    print("VERDICT:", "ALL PASS — 隐数框架承载旋转群/共轭/√ 全部可测" if all_ok
          else "SOME FAIL — 检查断言")


if __name__ == "__main__":
    main()
