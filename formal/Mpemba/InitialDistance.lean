import Mpemba.Relaxation

noncomputable section
namespace Mpemba

theorem hot_fourth_probability_lower :
    (248/1000 : ℝ) < gibbs (1/1000) 3 := by
  obtain ⟨h1l, h1u⟩ := hot_root_bounds 1
  obtain ⟨h2l, h2u⟩ := hot_root_bounds 2
  obtain ⟨h3l, h3u⟩ := hot_root_bounds 3
  norm_num [hotLower, hotUpper, q, Matrix.cons_val_two, Matrix.cons_val_three,
    Matrix.vecHead, Matrix.vecTail] at h1l h1u h2l h2u h3l h3u
  apply (lt_div_iff₀ (partition_positive (1/1000))).mpr
  norm_num [partition, q, Fin.sum_univ_succ, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.vecHead, Matrix.vecTail]
  linarith

def coldProbabilityUpper : State → ℝ := ![26/100, 26/100, 26/100, 24/100]

theorem cold_probability_upper (i : State) :
    gibbs (1/100) i ≤ coldProbabilityUpper i := by
  obtain ⟨h1l, h1u⟩ := cold_root_bounds 1
  obtain ⟨h2l, h2u⟩ := cold_root_bounds 2
  obtain ⟨h3l, h3u⟩ := cold_root_bounds 3
  norm_num [coldLower, coldUpper, q, Matrix.cons_val_two, Matrix.cons_val_three,
    Matrix.vecHead, Matrix.vecTail] at h1l h1u h2l h2u h3l h3u
  unfold gibbs
  apply (div_le_iff₀ (partition_positive (1/100))).mpr
  fin_cases i <;>
    norm_num [partition, q, coldProbabilityUpper, Fin.sum_univ_succ,
      Matrix.cons_val_two, Matrix.cons_val_three, Matrix.vecHead, Matrix.vecTail] <;>
    linarith

theorem bath_positive (i : State) : 0 < bath i := by
  unfold bath
  apply div_pos (q_positive i)
  norm_num [q, Fin.sum_univ_succ]

theorem cold_initial_distance_upper : bathDistance (gibbs (1/100)) < 24 := by
  have hbound (i : State) : (gibbs (1/100) i - bath i)^2 / q i ≤
      ((coldProbabilityUpper i)^2 + (bath i)^2) / q i := by
    apply div_le_div_of_nonneg_right _ (q_positive i).le
    have hp := (gibbs_positive (1/100) i).le
    have hu := cold_probability_upper i
    have hb := (bath_positive i).le
    have hs : (gibbs (1/100) i)^2 ≤ (coldProbabilityUpper i)^2 := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hu) (add_nonneg hp (le_trans hp hu))]
    nlinarith [mul_nonneg hp hb]
  have hsum := Finset.sum_le_sum (s := Finset.univ) (fun i _ => hbound i)
  change bathDistance (gibbs (1/100)) ≤ _ at hsum
  apply lt_of_le_of_lt hsum
  norm_num [coldProbabilityUpper, bath, q, Fin.sum_univ_succ]

theorem hot_initial_distance_lower : 24 < bathDistance (gibbs (1/1000)) := by
  have hterm : (gibbs (1/1000) 3 - bath 3)^2 / q 3 ≤
      bathDistance (gibbs (1/1000)) := by
    unfold bathDistance weightedSquare
    exact Finset.single_le_sum (f := fun i => (gibbs (1/1000) i - bath i)^2 / q i) (fun i _ =>
      div_nonneg (sq_nonneg _) (q_positive i).le) (Finset.mem_univ 3)
  have hl := hot_fourth_probability_lower
  have hlarge : 24 < (gibbs (1/1000) 3 - bath 3)^2 / q 3 := by
    have hq : q 3 = (1/400 : ℝ) := rfl
    have hb : bath 3 = (1/959 : ℝ) := by
      unfold bath
      rw [hq]
      norm_num [q, Fin.sum_univ_succ]
    rw [hq, hb]
    norm_num
    nlinarith [sq_nonneg (gibbs (1/1000) 3 - 248/1000)]
  exact lt_of_lt_of_le hlarge hterm

/-- Initial-to-eventual reversal for the actual master-equation solution. -/
theorem companion_distance_reversal :
    bathDistance (evolve (gibbs (1/100)) 0) <
      bathDistance (evolve (gibbs (1/1000)) 0) ∧
    ∀ᶠ t : ℝ in Filter.atTop,
      bathDistance (evolve (gibbs (1/1000)) t) <
      bathDistance (evolve (gibbs (1/100)) t) := by
  constructor
  · have hi (β : ℝ) : evolve (gibbs β) 0 = gibbs β := funext (evolve_initial _)
    rw [hi, hi]
    exact lt_trans cold_initial_distance_upper hot_initial_distance_lower
  · exact companion_eventual_distance_reversal

end Mpemba
