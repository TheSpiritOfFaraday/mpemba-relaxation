import Mpemba.ModeShape
import Mpemba.QuarterGibbs
import Mpemba.OriginalTriplets

noncomputable section
open Matrix
namespace Mpemba.Quarter

/-- A prefix contrast has strictly negative expectation below the bath
    inverse temperature, proved using exact real powers. -/
theorem power_contrast_negative {a b β : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hab : b < a) (hβ : β < 1) : a^β-(a/b)*b^β < 0 := by
  have hh := Real.rpow_lt_rpow_of_exponent_lt ((one_lt_div hb).mpr hab) hβ
  rw [Real.rpow_one, Real.div_rpow ha.le hb.le] at hh
  have h := (div_lt_iff₀ (Real.rpow_pos_of_pos hb β)).mp hh
  linarith

def prefixNumerator (q : FourWeights) (β : ℝ) : Fin 3 → ℝ :=
  ![1-(1/q.t)*q.t^β,
    1+q.t^β-(q.Q2/q.u)*q.u^β,
    1+q.t^β+q.u^β-(q.Q3/q.v)*q.v^β]

theorem prefixNumerator_negative (q : FourWeights) (hq : q.Valid) (β : ℝ)
    (hβ : β < 1) (i : Fin 3) : prefixNumerator q β i < 0 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have ht1 := hq.2.2.2
  have hut := hq.2.2.1
  have hvu := hq.2.1
  have h₁ := power_contrast_negative (by norm_num : (0:ℝ)<1) ht ht1 hβ
  have h₂ := power_contrast_negative (by norm_num : (0:ℝ)<1) hu (lt_trans hut ht1) hβ
  have h₃ := power_contrast_negative ht hu hut hβ
  have h₄ := power_contrast_negative (by norm_num : (0:ℝ)<1) hv
    (lt_trans hvu (lt_trans hut ht1)) hβ
  have h₅ := power_contrast_negative ht hv (lt_trans hvu hut) hβ
  have h₆ := power_contrast_negative hu hv hvu hβ
  simp only [Real.one_rpow] at h₁ h₂ h₄
  fin_cases i
  · exact h₁
  · simp only [prefixNumerator, cons_val_zero', cons_val_succ', FourWeights.Q2, add_div]
    linarith
  · simp only [prefixNumerator, cons_val_zero', cons_val_succ', FourWeights.Q3, add_div]
    linarith

def modeOverlap (q : FourWeights) (x : Fin 3 → ℝ) (β : ℝ) : ℝ :=
  let w := helmert q*ᵥx
  fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) β

theorem modeOverlap_numerator (q : FourWeights) (x : Fin 3 → ℝ) (β : ℝ) :
    modeOverlap q x β =
      (∑ i, x i*prefixNumerator q β i)/(Z3 q.t q.u β+q.v^β) := by
  simp [modeOverlap, fullMean, helmert, Matrix.mulVec, dotProduct,
    Fin.sum_univ_three, prefixNumerator]
  congr 1
  ring

theorem nonnegative_mode_overlap_negative (q : FourWeights) (hq : q.Valid)
    (x : Fin 3 → ℝ) (hx : ∀ i, 0 ≤ x i) (hxn : x ≠ 0)
    (β : ℝ) (hβ : β < 1) : modeOverlap q x β < 0 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hsum : (∑ i, x i*prefixNumerator q β i) < 0 := by
    have hi (i : Fin 3) : x i*prefixNumerator q β i ≤ 0 :=
      mul_nonpos_of_nonneg_of_nonpos (hx i) (prefixNumerator_negative q hq β hβ i).le
    have hex : ∃ i, x i ≠ 0 := by
      by_contra hn
      apply hxn
      ext i
      by_contra hni
      exact hn ⟨i, hni⟩
    obtain ⟨i, hxi⟩ := hex
    have hstrict : x i*prefixNumerator q β i < 0 :=
      mul_neg_of_pos_of_neg (lt_of_le_of_ne (hx i) (Ne.symm hxi))
        (prefixNumerator_negative q hq β hβ i)
    have hlt := Finset.sum_lt_sum (s := Finset.univ)
      (f := fun j => x j*prefixNumerator q β j) (g := fun _ => (0:ℝ))
      (fun j _ => hi j) ⟨i, Finset.mem_univ i, hstrict⟩
    simpa only [Finset.sum_const_zero] using hlt
  rw [modeOverlap_numerator]
  exact div_neg_of_neg_of_pos hsum (add_pos (Z3_pos ht hu β) (Real.rpow_pos_of_pos hv β))

theorem nonnegative_mode_quarter_bound (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x) (hx : ∀ i, 0 ≤ x i)
    (hot cold : ℝ) (hh : 0 < hot) (hhc : hot ≤ cold) (hc : cold < 1) :
    |modeOverlap q x cold|/|modeOverlap q x hot| < 4/3 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨hw₁, hw₂, hw₃⟩ := nonnegative_mode_shape q hq r hr hcone lam x h hx
  exact quarter_bound_from_mode_shape ht hq.2.2.2.le hu hq.2.2.1.le hv hq.2.1
    hh hhc hw₁ hw₂ hw₃
    (nonnegative_mode_overlap_negative q hq x hx h.nonzero hot (lt_of_le_of_lt hhc hc))
    (nonnegative_mode_overlap_negative q hq x hx h.nonzero cold hc)

theorem modeOverlap_smul (q : FourWeights) (x : Fin 3 → ℝ) (k β : ℝ) :
    modeOverlap q (k • x) β = k*modeOverlap q x β := by
  simp only [modeOverlap_numerator, Pi.smul_apply, smul_eq_mul,
    mul_assoc, ← Finset.mul_sum, mul_div_assoc]

theorem simple_top_mode_overlap_nonzero (q : FourWeights) (hq : q.Valid)
    (r : FourRates) (hcone : TripletCone q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x) (β : ℝ) (hβ : β < 1) :
    modeOverlap q x β ≠ 0 := by
  have hoff := contrast_offdiagonal_nonnegative q hq r hcone
  have hn := nonnegative_mode_overlap_negative q hq (fun i => |x i|)
    (fun i => abs_nonneg _) h.abs_nonzero β hβ
  obtain ⟨k, hk, heq⟩ := h.abs_is_nonzero_multiple hoff
  intro hz
  rw [heq, modeOverlap_smul, hz, mul_zero] at hn
  exact (lt_irrefl 0) hn

/-- The upper bound is independent of the sign and normalization of the
    eigenmode. The PSD certificate states the usual variational top-mode
    condition; no mode-shape or overlap-sign hypothesis is assumed. -/
theorem simple_top_mode_quarter_bound (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (hot cold : ℝ) (hh : 0 < hot) (hhc : hot ≤ cold) (hc : cold < 1) :
    |modeOverlap q x cold|/|modeOverlap q x hot| < 4/3 := by
  have hoff := contrast_offdiagonal_nonnegative q hq r hcone
  have habs := nonnegative_mode_quarter_bound q hq r hr hcone lam (fun i => |x i|)
    (h.abs_mode hoff) (fun i => abs_nonneg _) hot cold hh hhc hc
  obtain ⟨k, hk, heq⟩ := h.abs_is_nonzero_multiple hoff
  rw [heq, modeOverlap_smul, modeOverlap_smul, abs_mul, abs_mul,
    mul_div_mul_left _ _ (abs_ne_zero.mpr hk)] at habs
  exact habs

theorem original_triplets_quarter_bound (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Positive) (htrip : allOriginalTripletsFail q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (hot cold : ℝ) (hh : 0 < hot) (hhc : hot ≤ cold) (hc : cold < 1) :
    |modeOverlap q x cold|/|modeOverlap q x hot| < 4/3 :=
  simple_top_mode_quarter_bound q hq r hr.nonnegative
    ((allOriginalTripletsFail_iff_cone q hq r hr).mp htrip) lam x h hot cold hh hhc hc

end Mpemba.Quarter
