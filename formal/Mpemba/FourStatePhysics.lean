import Mpemba.OriginalTriplets
import Mpemba.QuarterGibbs
import Mathlib.Analysis.SpecialFunctions.Pow.Real

noncomputable section
open Matrix
namespace Mpemba.Quarter

theorem bath_vector_positive (q : FourWeights) (hq : q.Valid) (i : Fin 4) : 0 < q.vector i := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  fin_cases i <;> simp only [FourWeights.vector, cons_val_zero', cons_val_succ'] <;> positivity

theorem four_generator_conservative (q : FourWeights) (r : FourRates) (j : Fin 4) :
    (∑ i, generator q r i j) = 0 := by
  fin_cases j <;> simp [generator, Fin.sum_univ_four] <;> ring

theorem four_generator_detailed_balance (q : FourWeights) (hq : q.Valid)
    (r : FourRates) (i j : Fin 4) :
    generator q r i j*q.vector j = generator q r j i*q.vector i := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  fin_cases i <;> fin_cases j <;> simp [generator, FourWeights.vector]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu]

theorem four_generator_positive_rates (q : FourWeights) (hq : q.Valid)
    (r : FourRates) (hr : r.Positive) (i j : Fin 4) (hij : i ≠ j) :
    0 < generator q r i j := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨ha, hb, hc, hd, he, hf⟩ := hr
  fin_cases i <;> fin_cases j <;> norm_num at hij <;> simp [generator] <;> positivity

theorem four_gibbs_physical (q : FourWeights) (hq : q.Valid) (β : ℝ) (i : Fin 4) :
    Real.exp (-β*(-Real.log (q.vector i))) = (q.vector i)^β := by
  rw [Real.rpow_def_of_pos (bath_vector_positive q hq i)]
  congr 1
  ring

theorem fullMean_physical (q : FourWeights) (hq : q.Valid) (w : Fin 4 → ℝ) (β : ℝ) :
    fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) β =
      (∑ i, w i*Real.exp (-β*(-Real.log (q.vector i)))) /
        (∑ i, Real.exp (-β*(-Real.log (q.vector i)))) := by
  simp_rw [four_gibbs_physical q hq β]
  simp [fullMean, Z3, FourWeights.vector, Fin.sum_univ_four]

end Mpemba.Quarter
