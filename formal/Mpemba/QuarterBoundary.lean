import Mpemba.QuarterBound
import Mathlib.Tactic.Positivity

/-! The eigenvalue-free boundary sign calculation in the first-three
arithmetic restriction. `ModeShape.lean` connects this polynomial to the
contrast matrix and spectral argument. -/
namespace Mpemba.Quarter

def boundaryN (q₁ q₂ q₃ a b d : ℝ) : ℝ :=
  2*a*q₂*q₃*(q₁-q₂) + d*q₁*q₃*(2*q₂+q₃) - b*q₂*q₃*(2*q₁+q₃)

def boundaryJ (q₁ q₂ q₃ c e : ℝ) : ℝ :=
  e*q₁*(2*q₂+q₃) - c*q₂*(2*q₁+q₃)

def boundaryT (q₁ q₂ q₃ q₄ a b c d e f : ℝ) : ℝ :=
  (c+e+f)*boundaryN q₁ q₂ q₃ a b d +
    q₄*(2*c*e*q₃*(q₁-q₂)+f*boundaryJ q₁ q₂ q₃ c e)

theorem boundaryJ_nonneg_of_c_le_e {q₁ q₂ q₃ c e : ℝ}
    (hq₁ : 0 ≤ q₁) (hq₃ : 0 ≤ q₃) (h12 : 0 ≤ q₁-q₂)
    (hq₂ : 0 ≤ q₂) (hc : 0 ≤ c) (hce : c ≤ e) :
    0 ≤ boundaryJ q₁ q₂ q₃ c e := by
  have heq : boundaryJ q₁ q₂ q₃ c e =
      (e-c)*q₁*(2*q₂+q₃) + c*q₃*(q₁-q₂) := by
    unfold boundaryJ
    ring
  rw [heq]
  have hdiff : 0 ≤ e-c := sub_nonneg.mpr hce
  positivity

theorem boundaryJ_lift {q₁ q₂ q₃ q₄ c e : ℝ}
    (hq₁ : 0 < q₁) (_hq₂ : 0 ≤ q₂) (hq₃ : 0 ≤ q₃) (hq₄ : 0 ≤ q₄)
    (h12 : 0 ≤ q₁-q₂) (h34 : 0 ≤ q₃-q₄) (he : 0 ≤ e)
    (hF : c*q₂*(2*q₁+q₄) ≤ e*q₁*(2*q₂+q₄)) :
    0 ≤ boundaryJ q₁ q₂ q₃ c e := by
  have hJ₄ : 0 ≤ e*q₁*(2*q₂+q₄)-c*q₂*(2*q₁+q₄) := sub_nonneg.mpr hF
  have hid : boundaryJ q₁ q₂ q₃ c e * (2*q₁+q₄) =
      (e*q₁*(2*q₂+q₄)-c*q₂*(2*q₁+q₄))*(2*q₁+q₃) +
      2*e*q₁*(q₁-q₂)*(q₃-q₄) := by
    unfold boundaryJ
    ring
  have hprod : 0 ≤ boundaryJ q₁ q₂ q₃ c e * (2*q₁+q₄) := by
    rw [hid]
    positivity
  exact nonneg_of_mul_nonneg_left hprod (by positivity : 0 < 2*q₁+q₄)

/-- All three reset-boundary faces, including zero-rate equality cases. -/
theorem boundaryT_nonnegative {q₁ q₂ q₃ q₄ a b c d e f : ℝ}
    (hq₁ : 0 < q₁) (hq₂ : 0 < q₂) (hq₃ : 0 < q₃) (hq₄ : 0 < q₄)
    (h12 : q₂ ≤ q₁) (h34 : q₄ ≤ q₃)
    (_ha : 0 ≤ a) (_hb : 0 ≤ b) (hc : 0 ≤ c) (hd : 0 ≤ d) (he : 0 ≤ e) (hf : 0 ≤ f)
    (hN : 0 ≤ boundaryN q₁ q₂ q₃ a b d)
    (hF : c*q₂*(2*q₁+q₄) ≤ e*q₁*(2*q₂+q₄)+2*a*q₂*(q₁-q₂))
    (hface : a = 0 ∨ b = 0 ∨ c = 0) :
    0 ≤ boundaryT q₁ q₂ q₃ q₄ a b c d e f := by
  have h12' : 0 ≤ q₁-q₂ := sub_nonneg.mpr h12
  have h34' : 0 ≤ q₃-q₄ := sub_nonneg.mpr h34
  have positive_if_J (hJ : 0 ≤ boundaryJ q₁ q₂ q₃ c e) :
      0 ≤ boundaryT q₁ q₂ q₃ q₄ a b c d e f := by
    unfold boundaryT
    positivity
  rcases hface with hza | hzb | hzc
  · apply positive_if_J
    apply boundaryJ_lift hq₁ hq₂.le hq₃.le hq₄.le h12' h34' he
    simpa only [hza, mul_zero, zero_mul, add_zero] using hF
  · by_cases hce : c ≤ e
    · exact positive_if_J (boundaryJ_nonneg_of_c_le_e hq₁.le hq₃.le h12' hq₂.le hc hce)
    · have hec : 0 ≤ c-e := (sub_pos.mpr (lt_of_not_ge hce)).le
      have hslack : 0 ≤ e*q₁*(2*q₂+q₄)+2*a*q₂*(q₁-q₂)-c*q₂*(2*q₁+q₄) :=
        sub_nonneg.mpr hF
      have hidentity : boundaryN q₁ q₂ q₃ a b d + q₄*boundaryJ q₁ q₂ q₃ c e =
          q₃*(e*q₁*(2*q₂+q₄)+2*a*q₂*(q₁-q₂)-c*q₂*(2*q₁+q₄)) +
          2*q₁*q₂*(q₃-q₄)*(c-e)+d*q₁*q₃*(2*q₂+q₃) := by
        unfold boundaryN boundaryJ
        rw [hzb]
        ring
      have hNJ : 0 ≤ boundaryN q₁ q₂ q₃ a b d + q₄*boundaryJ q₁ q₂ q₃ c e := by
        rw [hidentity]
        positivity
      have hT : boundaryT q₁ q₂ q₃ q₄ a b c d e f =
          (c+e)*boundaryN q₁ q₂ q₃ a b d +
          f*(boundaryN q₁ q₂ q₃ a b d+q₄*boundaryJ q₁ q₂ q₃ c e) +
          2*c*e*q₃*q₄*(q₁-q₂) := by
        unfold boundaryT
        ring
      rw [hT]
      positivity
  · apply positive_if_J
    unfold boundaryJ
    rw [hzc]
    simp only [zero_mul, sub_zero]
    positivity

end Mpemba.Quarter
