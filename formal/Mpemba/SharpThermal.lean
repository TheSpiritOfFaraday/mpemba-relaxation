import Mpemba.SharpFamily
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity

noncomputable section
open Matrix Filter Topology
namespace Mpemba.Quarter.Sharp

/-- L parametrizes finite positive preparations. Choosing epsilon=exp(-L²)
    makes the highest-state Gibbs weights exp(-2/L) and exp(-2L). -/
def parameter (L : ℝ) : ℝ := Real.exp (-(L^2))
def hot (L : ℝ) : ℝ := (L⁻¹)^3
def cold (L : ℝ) : ℝ := L⁻¹
def thermalOverlap (L β : ℝ) : ℝ :=
  fullMean (weights (parameter L)).t (weights (parameter L)).u (weights (parameter L)).v
    (mode (parameter L) 0) (mode (parameter L) 1) (mode (parameter L) 2) (mode (parameter L) 3) β

theorem parameter_pos (L : ℝ) : 0 < parameter L := Real.exp_pos _

theorem parameter_tendsto_zero : Tendsto parameter atTop (𝓝 (0:ℝ)) :=
  Real.tendsto_exp_neg_atTop_nhds_zero.comp (tendsto_pow_atTop (by decide : (2:ℕ) ≠ 0))

theorem hot_tendsto_zero : Tendsto hot atTop (𝓝 (0:ℝ)) := by
  simpa only [hot, zero_pow (by decide : (3:ℕ) ≠ 0)] using
    (tendsto_inv_atTop_zero : Tendsto (fun L : ℝ => L⁻¹) atTop (𝓝 0)).pow 3

theorem cold_tendsto_zero : Tendsto cold atTop (𝓝 (0:ℝ)) := tendsto_inv_atTop_zero

theorem temperatures_ordered (L : ℝ) (hL : 1 < L) :
    0 < hot L ∧ hot L < cold L ∧ cold L < 1 := by
  have hpos : 0 < L := by linarith
  have hinv : 0 < L⁻¹ := inv_pos.mpr hpos
  have hinv1 : L⁻¹ < 1 := (inv_lt_one₀ hpos).mpr hL
  have hsq : (L⁻¹)^2 < 1 := by nlinarith
  have hcube := mul_lt_mul_of_pos_right hsq hinv
  exact ⟨pow_pos hinv 3, by dsimp only [hot, cold]; nlinarith only [hcube], hinv1⟩

theorem tail_hot_identity (L : ℝ) (hL : L ≠ 0) :
    (weights (parameter L)).v ^ hot L = Real.exp (-2*L⁻¹) := by
  change ((Real.exp (-(L^2)))^2) ^ ((L⁻¹)^3) = _
  rw [← Real.exp_nat_mul, Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
  congr 1
  norm_num
  field_simp [hL]
  ring

theorem tail_cold_identity (L : ℝ) :
    (weights (parameter L)).v ^ cold L = Real.exp (-2*L) := by
  change ((Real.exp (-(L^2)))^2) ^ L⁻¹ = _
  rw [← Real.exp_nat_mul, Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
  congr 1
  norm_num
  by_cases hL : L = 0
  · simp [hL]
  · field_simp [hL]
    ring

theorem tail_hot_limit :
    Tendsto (fun L => (weights (parameter L)).v ^ hot L) atTop (𝓝 (1:ℝ)) := by
  have hh : Tendsto (fun L : ℝ => Real.exp (-2*L⁻¹)) atTop (𝓝 (1:ℝ)) := by
    have harg := (tendsto_const_nhds : Tendsto (fun _ : ℝ => (-2:ℝ)) atTop (𝓝 (-2))).mul
      (tendsto_inv_atTop_zero : Tendsto (fun L : ℝ => L⁻¹) atTop (𝓝 0))
    simpa using (Real.continuous_exp.tendsto (-2*0)).comp harg
  apply hh.congr'
  filter_upwards [eventually_gt_atTop (0:ℝ)] with L hL
  exact (tail_hot_identity L (ne_of_gt hL)).symm

theorem tail_cold_limit :
    Tendsto (fun L => (weights (parameter L)).v ^ cold L) atTop (𝓝 (0:ℝ)) := by
  have hh := Real.tendsto_exp_neg_atTop_nhds_zero.comp
    (Filter.Tendsto.const_mul_atTop (by norm_num : (0:ℝ)<2) tendsto_id)
  simpa only [tail_cold_identity, neg_mul] using hh

theorem first_weight_limit {β : ℝ → ℝ} (hβ : Tendsto β atTop (𝓝 0)) :
    Tendsto (fun L => (weights (parameter L)).t ^ β L) atTop (𝓝 (1:ℝ)) := by
  have ht : Tendsto (fun L => (weights (parameter L)).t) atTop (𝓝 (1:ℝ)) := by
    simpa only [weights, zero_div, sub_zero] using
      (tendsto_const_nhds : Tendsto (fun _ : ℝ => (1:ℝ)) atTop (𝓝 1)).sub
        (parameter_tendsto_zero.div_const 3)
  simpa only [Real.rpow_zero] using ht.rpow hβ (Or.inl (by norm_num : (1:ℝ) ≠ 0))

theorem third_weight_limit {β : ℝ → ℝ} (hβ : Tendsto β atTop (𝓝 0)) :
    Tendsto (fun L => (weights (parameter L)).u ^ β L) atTop (𝓝 (1:ℝ)) := by
  simpa only [weights, Real.rpow_zero] using
    (tendsto_const_nhds : Tendsto (fun _ : ℝ => (3/5:ℝ)) atTop (𝓝 (3/5))).rpow hβ
      (Or.inl (by norm_num : (3/5:ℝ) ≠ 0))

theorem middle_mode_limit :
    Tendsto (fun L => mode (parameter L) 2) atTop (𝓝 (-1/3:ℝ)) := by
  have hc : ContinuousAt (fun ε : ℝ => mode ε 2) 0 := by simp [mode]; fun_prop
  simpa [mode] using hc.tendsto.comp parameter_tendsto_zero

theorem thermal_overlap_limit {β : ℝ → ℝ} (hβ : Tendsto β atTop (𝓝 0))
    (tail : ℝ) (hden : 3+tail ≠ 0)
    (hlast : Tendsto (fun L => (weights (parameter L)).v ^ β L) atTop (𝓝 tail)) :
    Tendsto (fun L => thermalOverlap L (β L)) atTop (𝓝 ((-2/15)/(3+tail))) := by
  have ht := first_weight_limit hβ
  have hu := third_weight_limit hβ
  have hn := (((tendsto_const_nhds : Tendsto (fun _ : ℝ => (1:ℝ)) atTop (𝓝 1)).add
    ((tendsto_const_nhds : Tendsto (fun _ : ℝ => (-4/5:ℝ)) atTop (𝓝 (-4/5))).mul ht)).add
    (middle_mode_limit.mul hu)).add (parameter_tendsto_zero.neg.mul hlast)
  have hd := (((tendsto_const_nhds : Tendsto (fun _ : ℝ => (1:ℝ)) atTop (𝓝 1)).add ht).add hu).add hlast
  have hdiv := hn.div hd (by norm_num; exact hden)
  convert hdiv using 1
  simp [thermalOverlap, fullMean, Z3, mode]
  ring

theorem hot_overlap_limit : Tendsto (fun L => thermalOverlap L (hot L)) atTop (𝓝 (-1/30:ℝ)) := by
  have hh := thermal_overlap_limit hot_tendsto_zero 1 (by norm_num) tail_hot_limit
  norm_num at hh
  simpa only [neg_div] using hh

theorem cold_overlap_limit : Tendsto (fun L => thermalOverlap L (cold L)) atTop (𝓝 (-2/45:ℝ)) := by
  have hh := thermal_overlap_limit cold_tendsto_zero 0 (by norm_num) tail_cold_limit
  norm_num at hh
  simpa only [neg_div] using hh

def advantage (L : ℝ) : ℝ := 1-|thermalOverlap L (hot L)|/|thermalOverlap L (cold L)|

theorem advantage_limit : Tendsto advantage atTop (𝓝 (1/4:ℝ)) := by
  have hh := (tendsto_const_nhds : Tendsto (fun _ : ℝ => (1:ℝ)) atTop (𝓝 1)).sub
    (hot_overlap_limit.abs.div cold_overlap_limit.abs (by norm_num))
  norm_num at hh
  have hn : (1:ℝ)-|(1/30:ℝ)|/|(2/45:ℝ)| = 1/4 := by
    rw [abs_of_pos (by norm_num : (0:ℝ)<1/30), abs_of_pos (by norm_num : (0:ℝ)<2/45)]
    norm_num
  rw [hn] at hh
  exact hh

end Mpemba.Quarter.Sharp
