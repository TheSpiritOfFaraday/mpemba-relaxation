import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! The final algebraic step of the quarter bound. `SlowModeBridge.lean`
derives the full-generator bound through the structural and Gibbs modules.
`SharpQuarter.lean` proves sharpness by a valid limiting family. -/
namespace Mpemba

/-- Given the conditional mean ordering and the first-three mass bounds,
    the four-to-three normalization factor bounds the amplitude ratio. -/
theorem quarter_bound_from_conditional_means
    (Ph Pc mh mc tail bh bc : ℝ)
    (hPh : 3/4 < Ph) (hPc : Pc < 1) (hP : Ph ≤ Pc)
    (hm : mh ≤ mc) (ht : tail ≤ 0)
    (hh : bh = Ph*mh + (1-Ph)*tail)
    (hc : bc = Pc*mc + (1-Pc)*tail)
    (hbh : bh < 0) (hbc : bc < 0) : |bc| / |bh| < 4/3 := by
  have hPh0 : 0 < Ph := by linarith
  have hPc0 : 0 < Pc := lt_of_lt_of_le hPh0 hP
  have hmeans : 0 ≤ Ph*Pc*(mc-mh) :=
    mul_nonneg (mul_pos hPh0 hPc0).le (sub_nonneg.mpr hm)
  have htail : 0 ≤ (Ph-Pc)*tail := mul_nonneg_of_nonpos_of_nonpos
    (sub_nonpos.mpr hP) ht
  have hcompare : Pc*bh ≤ Ph*bc := by
    rw [hh, hc]
    nlinarith only [hmeans, htail]
  have hfactor : Pc < (4/3)*Ph := by linarith
  have hstrict : (4/3)*Ph*bh < Pc*bh :=
    mul_lt_mul_of_neg_right hfactor hbh
  rw [abs_of_neg hbc, abs_of_neg hbh]
  apply (div_lt_iff₀ (neg_pos.mpr hbh)).mpr
  nlinarith only [hcompare, hstrict, hPh0]

theorem amplitude_advantage_below_quarter
    (ah ac : ℝ) (hh : 0 < ah) (hc : 0 < ac) (hr : ac / ah < 4/3) :
    1 - ah / ac < 1/4 := by
  have h := (div_lt_iff₀ hh).mp hr
  have h' : (3/4)*ac < ah := by linarith
  have h'' := (lt_div_iff₀ hc).mpr h'
  linarith

end Mpemba
