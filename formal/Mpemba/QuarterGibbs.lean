import Mpemba.Thermal
import Mpemba.QuarterBound
import Mathlib.Tactic.FieldSimp

noncomputable section
namespace Mpemba.Quarter

def Z3 (a b β : ℝ) : ℝ := 1 + a^β + b^β
def contrast₁ (a b β : ℝ) : ℝ := (1-a^β) / Z3 a b β
def contrast₂ (a b β : ℝ) : ℝ := (1+a^β-2*b^β) / Z3 a b β
def mean (a b w₀ w₁ w₂ β : ℝ) : ℝ := (w₀+w₁*a^β+w₂*b^β) / Z3 a b β

theorem Z3_pos {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (β : ℝ) : 0 < Z3 a b β := by
  unfold Z3
  linarith [Real.rpow_pos_of_pos ha β, Real.rpow_pos_of_pos hb β]

theorem contrast₁_monotone {a b h c : ℝ}
    (ha : 0 < a) (ha1 : a ≤ 1) (hb : 0 < b) (hb1 : b ≤ 1)
    (hh : 0 ≤ h) (hc : h ≤ c) : contrast₁ a b h ≤ contrast₁ a b c := by
  have hza := Real.rpow_le_rpow_of_exponent_ge ha ha1 hc
  have hzb := Real.rpow_le_rpow_of_exponent_ge hb hb1 hc
  have hz1 := Real.rpow_le_one ha.le ha1 (le_trans hh hc)
  unfold contrast₁
  apply div_le_div₀ (sub_nonneg.mpr hz1) (by linarith) (Z3_pos ha hb c)
  unfold Z3
  linarith

theorem third_probability_identity {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (β : ℝ) :
    b^β / Z3 a b β = 1 / ((1/b)^β+(a/b)^β+1) := by
  have hz := Z3_pos ha hb β
  have hs : 0 < (1/b)^β+(a/b)^β+1 := by
    linarith [Real.rpow_pos_of_pos (div_pos (by norm_num : (0 : ℝ) < 1) hb) β,
      Real.rpow_pos_of_pos (div_pos ha hb) β]
  apply (div_eq_div_iff (ne_of_gt hz) (ne_of_gt hs)).mpr
  rw [Real.div_rpow (by norm_num) hb.le, Real.div_rpow ha.le hb.le, Real.one_rpow]
  unfold Z3
  field_simp [ne_of_gt (Real.rpow_pos_of_pos hb β)]

theorem third_probability_decreases {a b h c : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hba : b ≤ a) (hb1 : b ≤ 1) (hc : h ≤ c) :
    b^c / Z3 a b c ≤ b^h / Z3 a b h := by
  rw [third_probability_identity ha hb c, third_probability_identity ha hb h]
  apply one_div_le_one_div_of_le
  · linarith [Real.rpow_pos_of_pos (div_pos (by norm_num : (0 : ℝ) < 1) hb) h,
      Real.rpow_pos_of_pos (div_pos ha hb) h]
  · have h₁ := Real.rpow_le_rpow_of_exponent_le ((one_le_div hb).mpr hb1) hc
    have h₂ := Real.rpow_le_rpow_of_exponent_le ((one_le_div hb).mpr hba) hc
    linarith

theorem contrast₂_identity {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (β : ℝ) :
    contrast₂ a b β = 1-3*(b^β/Z3 a b β) := by
  unfold contrast₂
  field_simp [ne_of_gt (Z3_pos ha hb β)]
  unfold Z3
  ring

theorem contrast₂_monotone {a b h c : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hba : b ≤ a) (hb1 : b ≤ 1) (hc : h ≤ c) :
    contrast₂ a b h ≤ contrast₂ a b c := by
  rw [contrast₂_identity ha hb h, contrast₂_identity ha hb c]
  linarith [third_probability_decreases ha hb hba hb1 hc]

theorem mean_contrast_decomposition {a b : ℝ} (ha : 0 < a) (hb : 0 < b)
    (w₀ w₁ w₂ β : ℝ) :
    mean a b w₀ w₁ w₂ β = (w₀+w₁+w₂)/3 +
      ((w₀-w₁)/2)*contrast₁ a b β + ((w₀+w₁-2*w₂)/6)*contrast₂ a b β := by
  unfold mean contrast₁ contrast₂
  field_simp [ne_of_gt (Z3_pos ha hb β)]
  unfold Z3
  ring

/-- The conditional Gibbs mean increases under exactly the first-three
    arithmetic inequalities used in the manuscript's upper-bound proof. -/
theorem conditional_mean_monotone {a b h c w₀ w₁ w₂ : ℝ}
    (ha : 0 < a) (ha1 : a ≤ 1) (hb : 0 < b) (hba : b ≤ a)
    (hh : 0 ≤ h) (hc : h ≤ c) (hw₁ : w₁ ≤ w₀) (hw₂ : 2*w₂ ≤ w₀+w₁) :
    mean a b w₀ w₁ w₂ h ≤ mean a b w₀ w₁ w₂ c := by
  have hb1 := le_trans hba ha1
  rw [mean_contrast_decomposition ha hb w₀ w₁ w₂ h,
    mean_contrast_decomposition ha hb w₀ w₁ w₂ c]
  have h₁ := mul_le_mul_of_nonneg_left (contrast₁_monotone ha ha1 hb hb1 hh hc)
    (by linarith : 0 ≤ (w₀-w₁)/2)
  have h₂ := mul_le_mul_of_nonneg_left (contrast₂_monotone ha hb hba hb1 hc)
    (by linarith : 0 ≤ (w₀+w₁-2*w₂)/6)
  linarith

def relativeSum (a b d β : ℝ) : ℝ := (1/d)^β+(a/d)^β+(b/d)^β
def mass (a b d β : ℝ) : ℝ := Z3 a b β / (Z3 a b β+d^β)
def fullMean (a b d w₀ w₁ w₂ w₃ β : ℝ) : ℝ :=
  (w₀+w₁*a^β+w₂*b^β+w₃*d^β) / (Z3 a b β+d^β)

theorem relativeSum_pos {a b d : ℝ} (ha : 0 < a) (hb : 0 < b) (hd : 0 < d)
    (β : ℝ) : 0 < relativeSum a b d β := by
  unfold relativeSum
  linarith [Real.rpow_pos_of_pos (div_pos (by norm_num : (0 : ℝ) < 1) hd) β,
    Real.rpow_pos_of_pos (div_pos ha hd) β, Real.rpow_pos_of_pos (div_pos hb hd) β]

theorem mass_relative_identity {a b d : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hd : 0 < d) (β : ℝ) :
    mass a b d β = relativeSum a b d β / (relativeSum a b d β+1) := by
  have hz : 0 < Z3 a b β+d^β := add_pos (Z3_pos ha hb β) (Real.rpow_pos_of_pos hd β)
  have hs : 0 < relativeSum a b d β+1 := by linarith [relativeSum_pos ha hb hd β]
  unfold mass
  apply (div_eq_div_iff (ne_of_gt hz) (ne_of_gt hs)).mpr
  unfold relativeSum
  rw [Real.div_rpow (by norm_num) hd.le, Real.div_rpow ha.le hd.le,
    Real.div_rpow hb.le hd.le, Real.one_rpow]
  unfold Z3
  field_simp [ne_of_gt (Real.rpow_pos_of_pos hd β)]

theorem mass_monotone {a b d h c : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hd : 0 < d)
    (hd1 : d ≤ 1) (hda : d ≤ a) (hdb : d ≤ b) (hc : h ≤ c) :
    mass a b d h ≤ mass a b d c := by
  rw [mass_relative_identity ha hb hd h, mass_relative_identity ha hb hd c]
  have hs : relativeSum a b d h ≤ relativeSum a b d c := by
    unfold relativeSum
    have h₁ := Real.rpow_le_rpow_of_exponent_le ((one_le_div hd).mpr hd1) hc
    have h₂ := Real.rpow_le_rpow_of_exponent_le ((one_le_div hd).mpr hda) hc
    have h₃ := Real.rpow_le_rpow_of_exponent_le ((one_le_div hd).mpr hdb) hc
    linarith
  have hh : 0 < relativeSum a b d h+1 := by linarith [relativeSum_pos ha hb hd h]
  have hc' : 0 < relativeSum a b d c+1 := by linarith [relativeSum_pos ha hb hd c]
  apply (div_le_div_iff₀ hh hc').mpr
  nlinarith only [hs]

theorem mass_bounds {a b d β : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hd : 0 < d)
    (hd1 : d < 1) (hda : d ≤ a) (hdb : d ≤ b) (hβ : 0 < β) :
    3/4 < mass a b d β ∧ mass a b d β < 1 := by
  rw [mass_relative_identity ha hb hd β]
  have h₁ : 1 < (1/d)^β := by
    have hh := Real.rpow_lt_rpow (by norm_num : (0 : ℝ) ≤ 1)
      ((one_lt_div hd).mpr hd1) hβ
    simpa only [Real.one_rpow] using hh
  have h₂ : 1 ≤ (a/d)^β := by
    have hh := Real.rpow_le_rpow_of_exponent_le ((one_le_div hd).mpr hda) hβ.le
    simpa only [Real.rpow_zero] using hh
  have h₃ : 1 ≤ (b/d)^β := by
    have hh := Real.rpow_le_rpow_of_exponent_le ((one_le_div hd).mpr hdb) hβ.le
    simpa only [Real.rpow_zero] using hh
  have hs : 3 < relativeSum a b d β := by unfold relativeSum; linarith
  have hden : 0 < relativeSum a b d β+1 := by linarith
  constructor
  · apply (lt_div_iff₀ hden).mpr
    linarith
  · apply (div_lt_one hden).mpr
    linarith

theorem fullMean_decomposition {a b d : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hd : 0 < d) (w₀ w₁ w₂ w₃ β : ℝ) :
    fullMean a b d w₀ w₁ w₂ w₃ β =
      mass a b d β * mean a b w₀ w₁ w₂ β + (1-mass a b d β)*w₃ := by
  unfold fullMean mass mean
  field_simp [ne_of_gt (Z3_pos ha hb β),
    ne_of_gt (add_pos (Z3_pos ha hb β) (Real.rpow_pos_of_pos hd β))]
  ring

/-- The analytic quarter bound from the mode-shape lemma. `ModeShape.lean`
    derives the shape hypotheses, and `CoolingOverlap.lean` the signs. -/
theorem quarter_bound_from_mode_shape {a b d h c w₀ w₁ w₂ w₃ : ℝ}
    (ha : 0 < a) (ha1 : a ≤ 1) (hb : 0 < b) (hba : b ≤ a)
    (hd : 0 < d) (hdb : d < b) (hh : 0 < h) (hc : h ≤ c)
    (hw₁ : w₁ ≤ w₀) (hw₂ : 2*w₂ ≤ w₀+w₁) (ht : w₃ ≤ 0)
    (hhot : fullMean a b d w₀ w₁ w₂ w₃ h < 0)
    (hcold : fullMean a b d w₀ w₁ w₂ w₃ c < 0) :
    |fullMean a b d w₀ w₁ w₂ w₃ c| / |fullMean a b d w₀ w₁ w₂ w₃ h| < 4/3 := by
  have hd1 : d < 1 := lt_of_lt_of_le hdb (le_trans hba ha1)
  have hda : d ≤ a := le_trans hdb.le hba
  apply Mpemba.quarter_bound_from_conditional_means
    (mass a b d h) (mass a b d c) (mean a b w₀ w₁ w₂ h) (mean a b w₀ w₁ w₂ c) w₃
  · exact (mass_bounds ha hb hd hd1 hda hdb.le hh).1
  · exact (mass_bounds ha hb hd hd1 hda hdb.le (lt_of_lt_of_le hh hc)).2
  · exact mass_monotone ha hb hd hd1.le hda hdb.le hc
  · exact conditional_mean_monotone ha ha1 hb hba hh.le hc hw₁ hw₂
  · exact ht
  · exact fullMean_decomposition ha hb hd w₀ w₁ w₂ w₃ h
  · exact fullMean_decomposition ha hb hd w₀ w₁ w₂ w₃ c
  · exact hhot
  · exact hcold

end Mpemba.Quarter
