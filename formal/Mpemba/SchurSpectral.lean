import Mpemba.TopMode
import Mpemba.QuarterSchur
import Mathlib.Tactic.FieldSimp

noncomputable section
open Matrix
namespace Mpemba.Quarter

def schurEntry (B : Matrix (Fin 3) (Fin 3) ℝ) (lam : ℝ) (i j : Fin 2) : ℝ :=
  B i.castSucc j.castSucc + B i.castSucc 2 * B 2 j.castSucc / (lam-B 2 2)

theorem simple_top_schur_diagonal_gap
    (η : Fin 3 → ℝ) (B : Matrix (Fin 3) (Fin 3) ℝ) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode η B lam x) (hx : 0 < x 0) :
    0 < lam-schurEntry B lam 1 1 := by
  have hδ := h.diagonal_gap 2 0 (by decide) (ne_of_gt hx)
  let y : Fin 3 → ℝ := ![0, 1, B 2 1/(lam-B 2 2)]
  have hid : y ⬝ᵥ weightedDefect η B lam *ᵥ y =
      η 1 * (lam-schurEntry B lam 1 1) := by
    simp [y, dotProduct, Matrix.mulVec, Fin.sum_univ_three, weightedDefect_apply, schurEntry]
    field_simp [ne_of_gt hδ]
    ring
  have hquad : 0 ≤ y ⬝ᵥ weightedDefect η B lam *ᵥ y := by
    simpa only [star_trivial, RCLike.re_to_real] using h.top.re_dotProduct_nonneg y
  rw [hid] at hquad
  have hnonneg : 0 ≤ lam-schurEntry B lam 1 1 := by
    by_contra hn
    have hh := mul_neg_of_pos_of_neg (h.weights_pos 1) (lt_of_not_ge hn)
    linarith
  apply lt_of_le_of_ne hnonneg
  intro heq
  have hzero : y ⬝ᵥ weightedDefect η B lam *ᵥ y = 0 := by
    rw [hid, ← heq, mul_zero]
  have hkernel : weightedDefect η B lam *ᵥ y = 0 :=
    (h.top.dotProduct_mulVec_zero_iff y).mp (by simpa only [star_trivial] using hzero)
  have heigen := (weightedDefect_kernel_iff η h.weights_pos B lam y).mp hkernel
  obtain ⟨c, hc⟩ := h.simple y heigen
  have hc0 : c = 0 := by
    have hi := congrFun hc 0
    change 0 = c*x 0 at hi
    exact (mul_eq_zero.mp hi.symm).resolve_right (ne_of_gt hx)
  have hi := congrFun hc 1
  change 1 = c*x 1 at hi
  rw [hc0, zero_mul] at hi
  norm_num at hi

theorem schur_eigen_equations
    (B : Matrix (Fin 3) (Fin 3) ℝ) (lam : ℝ) (x : Fin 3 → ℝ)
    (heigen : B *ᵥ x = lam • x) (hδ : lam-B 2 2 ≠ 0) :
    schurEntry B lam 0 0*x 0 + schurEntry B lam 0 1*x 1 = lam*x 0 ∧
    schurEntry B lam 1 0*x 0 + schurEntry B lam 1 1*x 1 = lam*x 1 := by
  have h₀ := congrFun heigen 0
  have h₁ := congrFun heigen 1
  have h₂ := congrFun heigen 2
  simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three, Pi.smul_apply,
    smul_eq_mul] at h₀ h₁ h₂
  have hx₂ : x 2 = (B 2 0*x 0+B 2 1*x 1)/(lam-B 2 2) := by
    apply (eq_div_iff hδ).mpr
    linarith only [h₂]
  constructor
  · calc
      _ = B 0 0*x 0+B 0 1*x 1+B 0 2*((B 2 0*x 0+B 2 1*x 1)/(lam-B 2 2)) := by
        simp only [schurEntry, Fin.castSucc_zero, Fin.castSucc_one]
        ring
      _ = lam*x 0 := by rwa [← hx₂]
  · calc
      _ = B 1 0*x 0+B 1 1*x 1+B 1 2*((B 2 0*x 0+B 2 1*x 1)/(lam-B 2 2)) := by
        simp only [schurEntry, Fin.castSucc_zero, Fin.castSucc_one]
        ring
      _ = lam*x 1 := by rwa [← hx₂]

theorem schur_ratio_comparison
    (η : Fin 3 → ℝ) (B : Matrix (Fin 3) (Fin 3) ℝ) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode η B lam x) (hoff : ∀ i j, i ≠ j → 0 ≤ B i j)
    (hx : 0 < x 0) (s : ℝ) (hs : 0 ≤ s)
    (hcone : 0 ≤ schurEntry B lam 1 0 + s*(schurEntry B lam 1 1-schurEntry B lam 0 0)-
      s^2*schurEntry B lam 0 1) : s*x 0 ≤ x 1 := by
  have hδ := h.diagonal_gap 2 0 (by decide) (ne_of_gt hx)
  obtain ⟨he₀, he₁⟩ := schur_eigen_equations B lam x h.eigen (ne_of_gt hδ)
  have h₀₁ : 0 ≤ schurEntry B lam 0 1 := by
    unfold schurEntry
    exact add_nonneg (hoff 0 1 (by decide))
      (div_nonneg (mul_nonneg (hoff 0 2 (by decide)) (hoff 2 1 (by decide))) hδ.le)
  have hstrict : 0 < lam-schurEntry B lam 1 1+schurEntry B lam 0 1*s :=
    add_pos_of_pos_of_nonneg (simple_top_schur_diagonal_gap η B lam x h hx)
      (mul_nonneg h₀₁ hs)
  exact eigenvector_ratio_comparison _ _ _ _ lam (x 0) (x 1) s hx he₀ he₁ hstrict hcone

end Mpemba.Quarter
