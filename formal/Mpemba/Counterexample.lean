import Mpemba.InitialDistance
import Mpemba.Multiplicity
import Mpemba.Uniqueness

noncomputable section
namespace Mpemba

/-- All characteristic roots relevant to real reversible relaxation lie no
    higher than the certified simple nonstationary mode. -/
theorem slow_mode_is_slowest (x : ℝ) (hx : M.charpoly.IsRoot x) (hz : x ≠ 0) :
    x ≤ slowRate := by
  by_cases hs : x = slowRate
  · exact hs.le
  · have hd : (Matrix.diagonal (fun _ : State => x) - M).det = 0 := by
      rwa [Polynomial.IsRoot, charpoly_eval_determinant] at hx
    have hf := other_characteristic_roots_faster x hd hz hs
    linarith [slowRate_location.1]

/-- A single theorem collecting the requested Appendix B proof obligations. -/
theorem appendixB_counterexample :
    (∀ i, 0 < q i) ∧
    (∀ i j, i < j → energy i < energy j) ∧
    (∀ i j, i ≠ j → 0 < M i j) ∧
    (∀ j, (∑ i, M i j) = 0) ∧
    (∀ i j, M i j * q j = M j i * q i) ∧
    (∀ i j k, i < j → j < k → ¬ originalTriplet i j k) ∧
    (∀ j, (∑ i, M i j * w i) = slowRate * w j) ∧
    M.charpoly.rootMultiplicity slowRate = 1 ∧
    (∀ x, M.charpoly.IsRoot x → x ≠ 0 → x ≤ slowRate) ∧
    (0 < (1/1000 : ℝ) ∧ (1/1000 : ℝ) < 1/100 ∧ (1/100 : ℝ) < 1) ∧
    |overlap (gibbs (1/1000))| < |overlap (gibbs (1/100))| ∧
    (bathDistance (evolve (gibbs (1/100)) 0) <
      bathDistance (evolve (gibbs (1/1000)) 0) ∧
    ∀ᶠ t : ℝ in Filter.atTop,
      bathDistance (evolve (gibbs (1/1000)) t) <
      bathDistance (evolve (gibbs (1/100)) t)) := by
  exact ⟨q_positive, energy_order, positive_rates, conservation, detailed_balance,
    all_original_triplets_fail, left_eigenvector, slow_algebraic_multiplicity_one,
    slow_mode_is_slowest, by norm_num, strict_thermal_overlap_reversal,
    companion_distance_reversal⟩

/-- The same reversal holds for any global solutions of the master equation
    with the stated Gibbs initial conditions, by uniqueness. -/
theorem master_equation_distance_reversal
    (hot cold : ℝ → State → ℝ)
    (hhot : ∀ t, HasDerivAt hot (act (hot t)) t)
    (hcold : ∀ t, HasDerivAt cold (act (cold t)) t)
    (hhot0 : hot 0 = gibbs (1/1000)) (hcold0 : cold 0 = gibbs (1/100)) :
    bathDistance (cold 0) < bathDistance (hot 0) ∧
    ∀ᶠ t : ℝ in Filter.atTop, bathDistance (hot t) < bathDistance (cold t) := by
  rw [evolve_unique _ hot hhot hhot0, evolve_unique _ cold hcold hcold0]
  exact companion_distance_reversal

theorem reversal_at_a_finite_positive_time : ∃ t : ℝ, 0 < t ∧
    bathDistance (evolve (gibbs (1/1000)) t) <
      bathDistance (evolve (gibbs (1/100)) t) := by
  exact ((Filter.eventually_gt_atTop (0 : ℝ)).and companion_eventual_distance_reversal).exists

end Mpemba
