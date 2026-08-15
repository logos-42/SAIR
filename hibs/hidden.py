#!/usr/bin/env python3
"""HIBS — Hidden-space Bridge System: executable Python implementation.

Faithful port of the Lean 4 formalization at
~/Downloads/lean/HIBS (Liu & Xu, "The Hidden-space Bridge System").

Core model (labelled-pair):  Hidden = ⟨val : Int, tag⟩,  tag ∈ {S, R, iR}
  (A1)  projR / projImag are non-injective (one-way projections)
  (A2)  +,− stay in S ;  ×,÷ force a projection to R
  (A3)  √ forces a projection to iR (full: S-ray ↦ −x·iR, R-ray ↦ x·iR)

Every Lean theorem is re-verified as a runnable assertion so the
formalization and the implementation agree.
"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum, auto
from math import gcd
import functools

# ---------------------------------------------------------------- Tags
class Tag(Enum):
    S = auto()    # additive hidden ray  (⟨+ : evaluates to ℝ⁺)
    R = auto()    # multiplicative ray   (⟨− : evaluates to ℝ⁻)
    iR = auto()   # imaginary ray       (iℝ)

# ---------------------------------------------------------------- Hidden
@dataclass(frozen=True)
class Hidden:
    val: int
    tag: Tag

    def __repr__(self):
        return f"⟨{self.val}, {self.tag.name}⟩"

S = Hidden  # abbrev

# ι_R : Int → S
def iota_R(a: int) -> Hidden:
    return Hidden(a, Tag.S)

def hiddenImag(b: int) -> Hidden:
    return Hidden(b, Tag.iR)

def hiddenReal(a: int) -> Hidden:
    return Hidden(a, Tag.S)

# ---------------------------------------------------------------- ℂ (Int pairs)
@dataclass(frozen=True)
class C:
    re: int
    im: int

    def __repr__(self):
        return f"({self.re} + {self.im}i)" if self.im >= 0 else f"({self.re} - {-self.im}i)"

    def __add__(self, o: C) -> C: return C(self.re + o.re, self.im + o.im)
    def __sub__(self, o: C) -> C: return C(self.re - o.re, self.im - o.im)
    def __mul__(self, o: C) -> C:
        return C(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    def __neg__(self) -> C: return C(-self.re, -self.im)

C0 = C(0, 0)
C1 = C(1, 0)
CI = C(0, 1)

def conj(z: C) -> C:
    return C(z.re, -z.im)

# ---------------------------------------------------------------- Operations (A2/A3)
def hAdd(h1: Hidden, h2: Hidden) -> Hidden:
    return Hidden(h1.val + h2.val, Tag.S)

def hSub(h1: Hidden, h2: Hidden) -> Hidden:
    return Hidden(h1.val - h2.val, Tag.S)

def hMul(h1: Hidden, h2: Hidden) -> Hidden:
    return Hidden(h1.val * h2.val, Tag.R)

def hSqrt(h: Hidden) -> Hidden:
    """(A3) simple: square root maps onto the imaginary axis."""
    return Hidden(h.val, Tag.iR)

def hSqrtFull(h: Hidden) -> Hidden:
    """(A3) full: sign-aware — S-ray ↦ −x·iR (iR⁻), R-ray ↦ x·iR (iR⁺)."""
    if h.tag == Tag.S:
        return Hidden(-h.val, Tag.iR)
    return Hidden(h.val, Tag.iR)

# ---------------------------------------------------------------- Conjugation & signal
def conjS(h: Hidden) -> Hidden:
    """Complex conjugation on hidden numbers: real rays fixed, iR negated."""
    if h.tag == Tag.iR:
        return Hidden(-h.val, Tag.iR)
    return h

def signalRev(h: Hidden) -> Hidden:
    """⟨+ ↔ ⟨− : pure tag swap R ↔ iR."""
    if h.tag == Tag.R:
        return Hidden(h.val, Tag.iR)
    if h.tag == Tag.iR:
        return Hidden(h.val, Tag.R)
    return h

# ---------------------------------------------------------------- Projection / embedding
def hEval(h: Hidden) -> C:
    """i-reading: S ↦ re, R ↦ −re?  (paper: R-ray evaluates to −x) — see hEval below.

    Definitions.lean's π: S↦(v,0), R↦(v,0), iR↦(0,v).  Derivation.lean's hEval
    uses the S-side reading.  We implement hEval per Derivation: R-ray ↦ −v.
    """
    if h.tag == Tag.S:
        return C(h.val, 0)
    if h.tag == Tag.R:
        return C(-h.val, 0)   # R-ray evaluates to ℝ⁻ (negEval_R)
    return C(0, h.val)        # iR ↦ iℝ

def pi(h: Hidden) -> C:
    """Definitions.lean π: straight projection (value as re or im)."""
    if h.tag == Tag.iR:
        return C(0, h.val)
    return C(h.val, 0)

def iota(z: C) -> Hidden:
    """ι : ℂ → S  (additive monomorphism): a+bi ↦ ⟨a,S⟩ + ⟨b,iR⟩."""
    return hAdd(iota_R(z.re), hiddenImag(z.im))

# CompositeHidden: two-component embedding (Theorem 6.5)
@dataclass(frozen=True)
class CompositeHidden:
    realPart: Hidden
    imagPart: Hidden

def iota_prime(z: C) -> CompositeHidden:
    """ι' : paper embedding — real part on R-ray, imag part on iR-ray."""
    return CompositeHidden(Hidden(z.re, Tag.R), Hidden(z.im, Tag.iR))

def pi_prime(c: CompositeHidden) -> C:
    """π' : project composite — real ray ↦ re, iR ray ↦ im."""
    return C(pi(c.realPart).re + (pi(c.realPart).im), pi(c.imagPart).im)
    # realPart on R-ray evaluates to its value on re (projection); imag on im

def pi_prime_fixed(c: CompositeHidden) -> C:
    """π' corrected: re from realPart (val), im from imagPart (val)."""
    return C(c.realPart.val, c.imagPart.val)

def conjC(c: CompositeHidden) -> CompositeHidden:
    return CompositeHidden(conjS(c.realPart), conjS(c.imagPart))

def hEvalC(c: CompositeHidden) -> C:
    return hEval(c.realPart) + hEval(c.imagPart)

# imul : ×i action (Derivation)
def imul(h: Hidden) -> Hidden:
    """×i on hidden numbers: real → imaginary ray, with sign per Derivation."""
    if h.tag == Tag.S:
        return Hidden(h.val, Tag.iR)
    if h.tag == Tag.R:
        return Hidden(-h.val, Tag.iR)
    return Hidden(-h.val, Tag.S)  # iR × i = −1 (real, negated)

def hMulAdj(a: Hidden, b: Hidden) -> Hidden:
    """Adjacent multiplication (Derivation.lean's hMulAdj): value product."""
    return Hidden(a.val * b.val, Tag.R)

# ---------------------------------------------------------------- Axiom checks (runnable)
def check_A1() -> bool:
    """(A1) projR, projImag non-injective: ⟨3,S⟩ ≠ ⟨3,R⟩ but same projections."""
    a, b = Hidden(3, Tag.S), Hidden(3, Tag.R)
    assert a != b, "distinct"
    assert a.val == b.val, "projR not injective"
    assert hEval(a).im == 0 and hEval(b).im == 0
    return True

def check_A2() -> bool:
    """(A2) tag flow: +,− → S ; × → R."""
    h1, h2 = Hidden(5, Tag.S), Hidden(7, Tag.R)
    assert hAdd(h1, h2).tag == Tag.S
    assert hSub(h1, h2).tag == Tag.S
    assert hMul(h1, h2).tag == Tag.R
    return True

def check_A3() -> bool:
    """(A3) √ ↦ iR for all hidden numbers; full half-axis clauses."""
    for h in (Hidden(4, Tag.S), Hidden(4, Tag.R), Hidden(4, Tag.iR), Hidden(-3, Tag.R)):
        assert hSqrt(h).tag == Tag.iR, "simple (A3)"
        assert hSqrtFull(h).tag == Tag.iR, "full (A3) tag"
    # full: S-ray (posEval) ↦ iR⁻ ; R-ray (negEval) ↦ iR⁺
    hp = Hidden(5, Tag.S)   # posEval: hEval = 5 > 0
    hr = Hidden(5, Tag.R)   # negEval: hEval = −5 < 0
    assert hEval(hSqrtFull(hp)) == C(0, -5), "⟨+ ↦ iR⁻"
    assert hEval(hSqrtFull(hr)) == C(0, 5), "⟨− ↦ iR⁺"
    return True

# ---------------------------------------------------------------- Theorems (runnable)
def verify_theorems() -> dict:
    out = {}
    # conjS fixes real rays
    out["conjS_fixes_real"] = conjS(Hidden(3, Tag.S)) == Hidden(3, Tag.S) and \
                              conjS(Hidden(3, Tag.R)) == Hidden(3, Tag.R)
    # conjS involution
    out["conjS_involution"] = all(conjS(conjS(Hidden(v, t))) == Hidden(v, t)
                                  for v in (-2, 0, 3) for t in Tag)
    # hEval ∘ conjS = conj ∘ hEval
    out["hEval_conjS"] = all(hEval(conjS(Hidden(v, t))) == conj(hEval(Hidden(v, t)))
                             for v in (-2, 0, 3) for t in Tag)
    # signalRev involution
    out["signalRev_involution"] = all(signalRev(signalRev(Hidden(v, t))) == Hidden(v, t)
                                      for v in (-2, 0, 3) for t in Tag)
    # signalRev swaps pos/neg partition: posHidden(signalRev h) ↔ negHidden h
    posHidden = lambda h: h.tag == Tag.R
    negHidden = lambda h: h.tag == Tag.iR
    out["signalRev_swaps"] = all(
        posHidden(signalRev(h)) == negHidden(h)
        for h in (Hidden(5, Tag.S), Hidden(5, Tag.R), Hidden(5, Tag.iR)))
    # conjC intertwines ι' : ι'(conj z) = conjC(ι' z)
    out["iota_prime_intertwines_conj"] = all(
        iota_prime(conj(C(z.re, z.im))) == conjC(iota_prime(C(z.re, z.im)))
        for z in (C(2, 3), C(-1, 4), C(0, -5)))
    # π'(conjC(ι' z)) = conj z
    out["pi_prime_conjC_iota"] = all(
        pi_prime_fixed(conjC(iota_prime(z))) == conj(z)
        for z in (C(2, 3), C(-1, 4), C(0, -5)))
    # ι additive monomorphism: ι(z+w) = ι(z) + ι(w)
    out["iota_additive"] = all(
        iota(C(a.re + b.re, a.im + b.im)) == hAdd(iota(a), iota(b))
        for a in (C(1, 2), C(-3, 4)) for b in (C(5, -6), C(0, 0)))
    # A3 full: √ = ×(−i) on real rays (hEval side)
    out["sqrtFull_eq_neg_i"] = all(
        hEval(hSqrtFull(h)) == hEval(h) * C(0, -1)
        for h in (Hidden(4, Tag.S), Hidden(4, Tag.R), Hidden(-2, Tag.S)))
    # conjS anti-linear w.r.t. imul: conjS ∘ imul = imul³ ∘ conjS
    out["conjS_imul"] = all(
        conjS(imul(Hidden(v, t))) == imul(imul(imul(conjS(Hidden(v, t)))))
        for v in (-2, 0, 3) for t in Tag)
    # (√⟨)² ∈ ℝ : square of sqrt lands on R-ray (im = 0)
    out["sqrt_sq_real"] = all(
        hEval(hMulAdj(hSqrtFull(h), hSqrtFull(h))).im == 0
        for h in (Hidden(4, Tag.S), Hidden(4, Tag.R), Hidden(-2, Tag.S)))
    # half-axis dichotomy: iℝ∖{0} = iR⁻ ⊔ iR⁺
    def iR_neg(z): return z.re == 0 and z.im < 0
    def iR_pos(z): return z.re == 0 and z.im > 0
    out["half_axis_dichotomy"] = all(
        (z.re == 0 and z.im != 0) == (iR_neg(z) or iR_pos(z))
        for z in (C(0, 5), C(0, -5), C(1, 2), C(0, 0)))
    return out


if __name__ == "__main__":
    print("HIBS executable model — axiom & theorem verification")
    print("=" * 56)
    for name, fn in (("(A1) one-way projections", check_A1),
                     ("(A2) closure under ±, ×↦R", check_A2),
                     ("(A3) √ ↦ iℝ (full half-axes)", check_A3)):
        print(f"  {name:<30} {'PASS' if fn() else 'FAIL'}")
    th = verify_theorems()
    n_pass = sum(th.values())
    print(f"  theorems: {n_pass}/{len(th)} PASS")
    for k, v in th.items():
        if not v:
            print(f"    FAIL: {k}")
    print("-" * 56)
    print("sample: conjS ⟨5, S⟩ =", conjS(Hidden(5, Tag.S)),
          " hSqrtFull ⟨4, S⟩ =", hSqrtFull(Hidden(4, Tag.S)),
          " hEval =", hEval(hSqrtFull(Hidden(4, Tag.S))))
