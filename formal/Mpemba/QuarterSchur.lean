import Mpemba.QuarterBound
import Mathlib.Tactic.FieldSimp

namespace Mpemba.Quarter

/-- A two-dimensional eigenvector ratio comparison without division by a
    possibly zero second coordinate or off-diagonal entry. -/
theorem eigenvector_ratio_comparison
    (a₁₁ a₁₂ a₂₁ a₂₂ lam x₁ x₂ s : ℝ)
    (hx₁ : 0 < x₁) (he₁ : a₁₁*x₁+a₁₂*x₂ = lam*x₁)
    (he₂ : a₂₁*x₁+a₂₂*x₂ = lam*x₂)
    (hstrict : 0 < lam-a₂₂+a₁₂*s)
    (hcone : 0 ≤ a₂₁+s*(a₂₂-a₁₁)-s^2*a₁₂) : s*x₁ ≤ x₂ := by
  have hid : (a₂₁+s*(a₂₂-a₁₁)-s^2*a₁₂)*x₁ =
      (x₂-s*x₁)*(lam-a₂₂+a₁₂*s) := by
    have hm := congrArg (fun z : ℝ => s*z) he₁
    dsimp only at hm
    nlinarith only [hm, he₂]
  have hprod : 0 ≤ (x₂-s*x₁)*(lam-a₂₂+a₁₂*s) := by
    rw [← hid]
    exact mul_nonneg hcone hx₁.le
  have hdiff : 0 ≤ x₂-s*x₁ := nonneg_of_mul_nonneg_left hprod hstrict
  linarith

/-- Exact Schur-complement identity used in the three-dimensional contrast proof. -/
theorem schur_ratio_numerator
    (b₁₁ b₁₂ b₁₃ b₂₁ b₂₂ b₂₃ b₃₁ b₃₂ δ s : ℝ) (hδ : δ ≠ 0) :
    δ*((b₂₁+b₂₃*b₃₁/δ)+s*((b₂₂+b₂₃*b₃₂/δ)-(b₁₁+b₁₃*b₃₁/δ))-
      s^2*(b₁₂+b₁₃*b₃₂/δ)) =
    δ*(b₂₁+s*(b₂₂-b₁₁)-s^2*b₁₂)+(b₂₃-s*b₁₃)*(b₃₁+s*b₃₂) := by
  field_simp [hδ]
  ring

theorem schur_sign_from_boundary (lam b₃₃ C D : ℝ)
    (hlam : lam ≤ 0) (hδ : 0 ≤ lam-b₃₃) (hD : 0 ≤ D)
    (hboundary : 0 ≤ -b₃₃*C+D) : 0 ≤ (lam-b₃₃)*C+D := by
  by_cases hC : 0 ≤ C
  · exact add_nonneg (mul_nonneg hδ hC) hD
  · have hlamC := mul_nonneg_of_nonpos_of_nonpos hlam (lt_of_not_ge hC).le
    nlinarith only [hlamC, hboundary]

end Mpemba.Quarter
