import Mpemba.Data
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Algebra.BigOperators.Field

noncomputable section
namespace Mpemba

def energy (i : State) : ℝ := -Real.log (q i)

def partition (β : ℝ) : ℝ := ∑ i, (q i) ^ β

def gibbs (β : ℝ) (i : State) : ℝ := (q i) ^ β / partition β

theorem gibbs_weight_physical (β : ℝ) (i : State) :
    Real.exp (-β * energy i) = (q i) ^ β := by
  rw [Real.rpow_def_of_pos (q_positive i)]
  congr 1
  unfold energy
  ring

theorem energy_order (i j : State) (h : i < j) : energy i < energy j := by
  exact neg_lt_neg (Real.log_lt_log (q_positive j) (q_strictly_decreasing i j h))

theorem partition_positive (β : ℝ) : 0 < partition β := by
  unfold partition
  exact Finset.sum_pos (fun i _ => Real.rpow_pos_of_pos (q_positive i) β)
    Finset.univ_nonempty

theorem gibbs_positive (β : ℝ) (i : State) : 0 < gibbs β i :=
  div_pos (Real.rpow_pos_of_pos (q_positive i) β) (partition_positive β)

theorem gibbs_sum (β : ℝ) : (∑ i, gibbs β i) = 1 := by
  simp only [gibbs, ← Finset.sum_div]
  exact div_self (ne_of_gt (partition_positive β))

theorem gibbs_bath (i : State) : gibbs 1 i = q i / (∑ j, q j) := by
  simp [gibbs, partition]

theorem gibbs_bath_overlap_zero : overlap (gibbs 1) = 0 := by
  simp_rw [overlap, gibbs_bath, ← mul_div_assoc]
  rw [← Finset.sum_div]
  have hm : (∑ i, w i * q i) = 0 := by
    simpa only [mul_comm] using zero_bath_mean
  rw [hm, zero_div]

theorem root_enclosure {r l u : ℝ} {n : ℕ} (hr : 0 ≤ r)
    (hl : 0 ≤ l) (hn : n ≠ 0) (hlo : l ^ n ≤ r) (hhi : r ≤ u ^ n)
    (hu : 0 ≤ u) : l ≤ r ^ (n : ℝ)⁻¹ ∧ r ^ (n : ℝ)⁻¹ ≤ u := by
  have hp := Real.rpow_nonneg hr ((n : ℝ)⁻¹)
  constructor
  · apply (pow_le_pow_iff_left₀ hl hp hn).mp
    simpa only [Real.rpow_inv_natCast_pow hr hn] using hlo
  · apply (pow_le_pow_iff_left₀ hp hu hn).mp
    simpa only [Real.rpow_inv_natCast_pow hr hn] using hhi

def hotLower : State → ℝ := ![1, 999994/1000000, 999084/1000000, 994026/1000000]
def hotUpper : State → ℝ := ![1, 999995/1000000, 999085/1000000, 994027/1000000]
def coldLower : State → ℝ := ![1, 999949/1000000, 990878/1000000, 941844/1000000]
def coldUpper : State → ℝ := ![1, 999950/1000000, 990879/1000000, 941845/1000000]

set_option maxRecDepth 10000 in
set_option maxHeartbeats 0 in
theorem hot_root_bounds (i : State) :
    hotLower i ≤ (q i) ^ (1/1000 : ℝ) ∧ (q i) ^ (1/1000 : ℝ) ≤ hotUpper i := by
  have hn : (1000 : ℕ) ≠ 0 := by norm_num
  have h := root_enclosure (n := 1000) (r := q i) (l := hotLower i) (u := hotUpper i)
    (q_positive i).le
  simp only [Nat.cast_ofNat, one_div] at h ⊢
  fin_cases i <;>
    apply h (by norm_num [hotLower]) hn (by norm_num [hotLower, q])
      (by norm_num [hotUpper, q]) (by norm_num [hotUpper])

set_option maxRecDepth 10000 in
set_option maxHeartbeats 0 in
theorem cold_root_bounds (i : State) :
    coldLower i ≤ (q i) ^ (1/100 : ℝ) ∧ (q i) ^ (1/100 : ℝ) ≤ coldUpper i := by
  have hn : (100 : ℕ) ≠ 0 := by norm_num
  have h := root_enclosure (n := 100) (r := q i) (l := coldLower i) (u := coldUpper i)
    (q_positive i).le
  simp only [Nat.cast_ofNat, one_div] at h ⊢
  fin_cases i <;>
    apply h (by norm_num [coldLower]) hn (by norm_num [coldLower, q])
      (by norm_num [coldUpper, q]) (by norm_num [coldUpper])

theorem overlap_formula (β : ℝ) : overlap (gibbs β) =
    (1 - (759/1000) * (q 1) ^ β - (48957/80000) * (q 2) ^ β -
      (1/250) * (q 3) ^ β) / (1 + (q 1) ^ β + (q 2) ^ β + (q 3) ^ β) := by
  simp [overlap, gibbs, partition, Fin.sum_univ_succ, w, q, sub_div, add_div]
  ring

theorem hot_overlap_bounds :
    -(9376/100000 : ℝ) < overlap (gibbs (1/1000)) ∧
    overlap (gibbs (1/1000)) < -(9375/100000 : ℝ) := by
  obtain ⟨h1l, h1u⟩ := hot_root_bounds 1
  obtain ⟨h2l, h2u⟩ := hot_root_bounds 2
  obtain ⟨h3l, h3u⟩ := hot_root_bounds 3
  norm_num [hotLower, hotUpper, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.vecHead, Matrix.vecTail] at h1l h1u h2l h2u h3l h3u
  have hd : 0 < 1 + (q 1) ^ (1/1000 : ℝ) + (q 2) ^ (1/1000 : ℝ) +
    (q 3) ^ (1/1000 : ℝ) := by linarith
  rw [overlap_formula]
  constructor
  · apply (lt_div_iff₀ hd).mpr
    linarith
  · apply (div_lt_iff₀ hd).mpr
    linarith

theorem cold_overlap_bounds :
    -(9386/100000 : ℝ) < overlap (gibbs (1/100)) ∧
    overlap (gibbs (1/100)) < -(9385/100000 : ℝ) := by
  obtain ⟨h1l, h1u⟩ := cold_root_bounds 1
  obtain ⟨h2l, h2u⟩ := cold_root_bounds 2
  obtain ⟨h3l, h3u⟩ := cold_root_bounds 3
  norm_num [coldLower, coldUpper, Matrix.cons_val_two, Matrix.cons_val_three, Matrix.vecHead, Matrix.vecTail] at h1l h1u h2l h2u h3l h3u
  have hd : 0 < 1 + (q 1) ^ (1/100 : ℝ) + (q 2) ^ (1/100 : ℝ) +
    (q 3) ^ (1/100 : ℝ) := by linarith
  rw [overlap_formula]
  constructor
  · apply (lt_div_iff₀ hd).mpr
    linarith
  · apply (div_lt_iff₀ hd).mpr
    linarith

theorem strict_thermal_overlap_reversal :
    |overlap (gibbs (1/1000))| < |overlap (gibbs (1/100))| := by
  obtain ⟨hl, hu⟩ := hot_overlap_bounds
  obtain ⟨cl, cu⟩ := cold_overlap_bounds
  rw [abs_of_neg (by linarith : overlap (gibbs (1/1000)) < 0),
    abs_of_neg (by linarith : overlap (gibbs (1/100)) < 0)]
  linarith

end Mpemba
