import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

namespace Mpemba

theorem branch_failure_cone (F G H h d : ℝ) (hh : 0 < h) (hd : 0 < d)
    (he : F = h*G + d*H) :
    ¬ ((0 ≤ H ∧ 0 < F) ∨ (H ≤ 0 ∧ 0 < G)) ↔ F ≤ 0 ∧ G ≤ 0 := by
  constructor
  · intro hn
    by_cases hH : 0 ≤ H
    · have hF : F ≤ 0 := le_of_not_gt (fun hF => hn (Or.inl ⟨hH, hF⟩))
      refine ⟨hF, le_of_not_gt ?_⟩
      intro hG
      nlinarith [mul_pos hh hG, mul_nonneg hd.le hH]
    · have hH' : H ≤ 0 := (lt_of_not_ge hH).le
      have hG : G ≤ 0 := le_of_not_gt (fun hG => hn (Or.inr ⟨hH', hG⟩))
      refine ⟨?_, hG⟩
      nlinarith [mul_nonpos_of_nonneg_of_nonpos hh.le hG,
        mul_nonpos_of_nonneg_of_nonpos hd.le hH']
  · rintro ⟨hF, hG⟩ (⟨_, hp⟩ | ⟨_, hp⟩)
    · exact (not_lt_of_ge hF) hp
    · exact (not_lt_of_ge hG) hp

/-- Full branch equivalence, retaining the inclusive equality case κ=1.
    This is the abstract form of the manuscript's first lemma. -/
theorem original_triplet_failure_iff (A B C t u : ℝ)
    (hC : 0 < C) (ht : 0 < t) (ht1 : t < 1) (hu : 0 < u) :
    let F := A*(1+t*u/2) - B*(1+u/2) - C*(1-t)
    let G := A*t-B
    let κ := (A+B)/(C*(1+t))
    ¬ ((1 ≤ κ ∧ 0 < F) ∨ (κ ≤ 1 ∧ 0 < G)) ↔ F ≤ 0 ∧ G ≤ 0 := by
  dsimp only
  have htplus : 0 < 1+t := by linarith
  have hden : 0 < C*(1+t) := mul_pos hC htplus
  have hk1 : 1 ≤ (A+B)/(C*(1+t)) ↔ 0 ≤ A+B-C*(1+t) := by
    rw [le_div_iff₀ hden]
    constructor <;> intro h <;> linarith
  have hk2 : (A+B)/(C*(1+t)) ≤ 1 ↔ A+B-C*(1+t) ≤ 0 := by
    rw [div_le_iff₀ hden]
    constructor <;> intro h <;> linarith
  rw [hk1, hk2]
  apply branch_failure_cone _ _ _ (2/(1+t)+u/2) ((1-t)/(1+t))
  · exact add_pos (div_pos (by norm_num) htplus) (div_pos hu (by norm_num))
  · exact div_pos (by linarith) htplus
  · field_simp [ne_of_gt htplus]
    ring

end Mpemba
