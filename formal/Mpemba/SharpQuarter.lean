import Mpemba.SharpThermal
import Mpemba.FourStatePhysics

noncomputable section
open Matrix Filter Topology
namespace Mpemba.Quarter

def amplitudeAdvantage (q : FourWeights) (w : Fin 4 → ℝ) (hot cold : ℝ) : ℝ :=
  1-|fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) hot| /
    |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) cold|

/-- Attainable ordinary cooling amplitude advantages at bath inverse
    temperature one, under all four original-rate triplet failures. -/
def attainableAdvantages : Set ℝ := {a | ∃ (q : FourWeights) (r : FourRates)
    (lam : ℝ) (w : Fin 4 → ℝ) (hot cold : ℝ),
    q.Valid ∧ r.Positive ∧ allOriginalTripletsFail q r ∧ SimpleSlowMode q r lam w ∧
    0 < hot ∧ hot < cold ∧ cold < 1 ∧
    |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) hot| <
      |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) cold| ∧
    a = amplitudeAdvantage q w hot cold}

theorem attainable_advantage_strict_upper_bound {a : ℝ} (ha : a ∈ attainableAdvantages) :
    a < 1/4 := by
  obtain ⟨q, r, lam, w, hot, cold, hq, hr, htrip, hslow, hh, hhc, hc, _, rfl⟩ := ha
  exact full_generator_amplitude_advantage_below_quarter q hq r hr htrip lam w hslow
    hot cold hh hhc.le hc

namespace Sharp

theorem eventually_amplitude_reversal :
    ∀ᶠ L in atTop, |thermalOverlap L (hot L)| < |thermalOverlap L (cold L)| := by
  have hdiff := cold_overlap_limit.abs.sub hot_overlap_limit.abs
  have hval : |(-2/45:ℝ)|-|(-1/30:ℝ)| = 1/90 := by
    rw [abs_of_neg (by norm_num : (-2/45:ℝ)<0), abs_of_neg (by norm_num : (-1/30:ℝ)<0)]
    norm_num
  rw [hval] at hdiff
  have hpos := (tendsto_order.mp hdiff).1 0 (by norm_num : (0:ℝ)<1/90)
  filter_upwards [hpos] with L hL
  linarith

/-- Every sufficiently large finite parameter gives a valid positive
    reversible generator, a simple slow mode and a strict cooling reversal. -/
theorem eventually_family_attainable : ∀ᶠ L in atTop, advantage L ∈ attainableAdvantages := by
  have hchecks := parameter_tendsto_zero.eventually eventually_checks_positive
  have hsmall := (tendsto_order.mp parameter_tendsto_zero).2 (1/30) (by norm_num : (0:ℝ)<1/30)
  filter_upwards [hchecks, hsmall, eventually_gt_atTop (1:ℝ), eventually_amplitude_reversal]
    with L hcheck hsize hL hrev
  have hq := weights_valid (parameter L) (parameter_pos L) hsize
  have hr := rates_positive_of_checks (parameter L) hcheck
  have hcone := cone_of_checks (parameter L) (parameter_pos L) hsize hcheck
  have htrip := (allOriginalTripletsFail_iff_cone _ hq _ hr).mpr hcone
  have hslow := family_simple_slow (parameter L) (parameter_pos L) hsize hcheck
  obtain ⟨hh, hhc, hc⟩ := temperatures_ordered L hL
  exact ⟨weights (parameter L), rates (parameter L), -3, mode (parameter L), hot L, cold L,
    hq, hr, htrip, hslow, hh, hhc, hc, hrev, rfl⟩

end Sharp

/-- The 25% number is the sharp supremum of actual admissible ordinary
    cooling amplitude advantages, rather than a conditional scalar bound. -/
theorem sharp_quarter_isLUB : IsLUB attainableAdvantages (1/4:ℝ) := by
  constructor
  · intro a ha
    exact (attainable_advantage_strict_upper_bound ha).le
  · intro b hb
    by_contra hnot
    have hbq : b < (1/4:ℝ) := lt_of_not_ge hnot
    have hnear := (tendsto_order.mp Sharp.advantage_limit).1 b hbq
    obtain ⟨L, hmem, hgt⟩ := (Sharp.eventually_family_attainable.and hnear).exists
    exact (not_le_of_gt hgt) (hb hmem)

theorem sharp_quarter_not_attained : (1/4:ℝ) ∉ attainableAdvantages := by
  intro h
  exact (lt_irrefl (1/4:ℝ)) (attainable_advantage_strict_upper_bound h)

theorem exists_advantage_arbitrarily_close (δ : ℝ) (hδ : 0 < δ) :
    ∃ a ∈ attainableAdvantages, (1/4:ℝ)-δ < a ∧ a < 1/4 := by
  have hnear := (tendsto_order.mp Sharp.advantage_limit).1 (1/4-δ) (by linarith)
  obtain ⟨L, hmem, hgt⟩ := (Sharp.eventually_family_attainable.and hnear).exists
  exact ⟨Sharp.advantage L, hmem, hgt, attainable_advantage_strict_upper_bound hmem⟩

end Mpemba.Quarter
