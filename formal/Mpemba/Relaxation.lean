import Mpemba.Dynamics

noncomputable section
open Filter Topology
namespace Mpemba

/-- A positive quadratic distance to the bath, proportional to Pearson chi-square. -/
def weightedSquare (x : State → ℝ) : ℝ := ∑ i, x i ^ 2 / q i

def bathDistance (p : State → ℝ) : ℝ := weightedSquare (fun i => p i - bath i)

def scaledDeviation (p : State → ℝ) (t : ℝ) (i : State) : ℝ :=
  overlap p * rightMode i +
    Real.exp ((fastRate₁ - slowRate) * t) * fastComponent₁ p i +
    Real.exp ((fastRate₂ - slowRate) * t) * fastComponent₂ p i

theorem evolve_deviation (p : State → ℝ) (hp : total p = 1) (t : ℝ) (i : State) :
    evolve p t i - bath i = Real.exp (slowRate*t) * scaledDeviation p t i := by
  have h₁ : Real.exp (slowRate*t) * Real.exp ((fastRate₁-slowRate)*t) =
      Real.exp (fastRate₁*t) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₂ : Real.exp (slowRate*t) * Real.exp ((fastRate₂-slowRate)*t) =
      Real.exp (fastRate₂*t) := by
    rw [← Real.exp_add]
    congr 1
    ring
  unfold evolve scaledDeviation
  rw [hp]
  simp only [mul_add, ← mul_assoc, h₁, h₂, one_mul]
  ring

theorem weighted_square_mul (a : ℝ) (x : State → ℝ) :
    weightedSquare (fun i => a * x i) = a^2 * weightedSquare x := by
  simp only [weightedSquare, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem distance_factorization (p : State → ℝ) (hp : total p = 1) (t : ℝ) :
    bathDistance (evolve p t) = Real.exp (slowRate*t)^2 *
      weightedSquare (scaledDeviation p t) := by
  unfold bathDistance
  simp_rw [evolve_deviation p hp t]
  exact weighted_square_mul _ _

theorem fast_relative_exp_limit {r : ℝ} (hr : r < slowRate) :
    Tendsto (fun t : ℝ => Real.exp ((r-slowRate)*t)) atTop (𝓝 0) := by
  exact Real.tendsto_exp_atBot.comp
    ((tendsto_const_mul_atBot_of_neg (sub_neg.mpr hr)).mpr tendsto_id)

theorem scaled_deviation_limit (p : State → ℝ) (i : State) :
    Tendsto (fun t => scaledDeviation p t i) atTop (𝓝 (overlap p * rightMode i)) := by
  have h₁ := fast_relative_exp_limit (r := fastRate₁)
    (lt_trans fast_rates_below.1 (by linarith [slowRate_location.1]))
  have h₂ := fast_relative_exp_limit (r := fastRate₂)
    (lt_trans fast_rates_below.2 (by linarith [slowRate_location.1]))
  have h := ((tendsto_const_nhds (x := overlap p * rightMode i)).add
    (h₁.mul_const (fastComponent₁ p i))).add (h₂.mul_const (fastComponent₂ p i))
  simpa only [zero_mul, add_zero, scaledDeviation] using h

theorem scaled_distance_limit (p : State → ℝ) :
    Tendsto (fun t => weightedSquare (scaledDeviation p t)) atTop
      (𝓝 ((overlap p)^2 * weightedSquare rightMode)) := by
  have h : Tendsto (fun t => ∑ i, (scaledDeviation p t i)^2 / q i) atTop
      (𝓝 (∑ i, (overlap p * rightMode i)^2 / q i)) := by
    apply tendsto_finset_sum
    intro i _
    exact ((scaled_deviation_limit p i).pow 2).div_const (q i)
  change Tendsto (fun t => weightedSquare (scaledDeviation p t)) atTop
    (𝓝 (weightedSquare (fun i => overlap p * rightMode i))) at h
  rwa [weighted_square_mul] at h

theorem right_mode_weighted_square_positive : 0 < weightedSquare rightMode := by
  norm_num [weightedSquare, rightMode, q, w, modeNorm, Fin.sum_univ_succ]

/-- Strict overlap ordering entails eventual ordering in a genuine bath distance
    for solutions of the printed generator, with no asymptotic premise assumed. -/
theorem overlap_implies_eventual_distance_order (p r : State → ℝ)
    (hp : total p = 1) (hr : total r = 1) (hamp : |overlap p| < |overlap r|) :
    ∀ᶠ t : ℝ in atTop, bathDistance (evolve p t) < bathDistance (evolve r t) := by
  have hsq : (overlap p)^2 < (overlap r)^2 := by
    nlinarith [sq_abs (overlap p), sq_abs (overlap r), abs_nonneg (overlap p)]
  have hlim : (overlap p)^2 * weightedSquare rightMode <
      (overlap r)^2 * weightedSquare rightMode :=
    mul_lt_mul_of_pos_right hsq right_mode_weighted_square_positive
  have hdiff := (scaled_distance_limit r).sub (scaled_distance_limit p)
  have he : ∀ᶠ t : ℝ in atTop,
      0 < weightedSquare (scaledDeviation r t) - weightedSquare (scaledDeviation p t) :=
    hdiff.eventually_const_lt (sub_pos.mpr hlim)
  filter_upwards [he] with t ht
  rw [distance_factorization p hp, distance_factorization r hr]
  exact mul_lt_mul_of_pos_left (sub_pos.mp ht) (sq_pos_of_pos (Real.exp_pos _))

theorem companion_eventual_distance_reversal :
    ∀ᶠ t : ℝ in atTop,
      bathDistance (evolve (gibbs (1/1000)) t) <
      bathDistance (evolve (gibbs (1/100)) t) := by
  apply overlap_implies_eventual_distance_order
  · exact gibbs_sum _
  · exact gibbs_sum _
  · exact strict_thermal_overlap_reversal

end Mpemba
