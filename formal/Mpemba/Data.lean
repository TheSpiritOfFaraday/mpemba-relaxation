import Mathlib.Data.Matrix.Notation
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

/-! Exact data printed in Appendix B. Indices 0,1,2,3 represent states 1,2,3,4.
    `M i j` is the rate from j to i (column probability convention). -/
namespace Mpemba

noncomputable section

abbrev State := Fin 4

def q : State → ℝ := ![1, 199/200, 2/5, 1/400]

def M : Matrix State State ℝ := !![
  -49029/4000, 1, 28, 249/10;
  199/200, -10689318350801/117041850000, 330154675693/1470375000, 4983/200;
  56/5, 660309351386/7315115625, -36129268114499561/143029257750000,
    1527327066/151990625;
  249/4000, 4983/79600, 763663533/12159250000, -72789920403/1215925000]

def w : State → ℝ := ![1, -759/1000, -48957/80000, -1/250]

def slowRate : ℝ := -4966671/250000

def modeNorm : ℝ := 27567998009/16000000000

def rightMode (i : State) : ℝ := q i * w i / modeNorm

def overlap (p : State → ℝ) : ℝ := ∑ i, w i * p i

def tripletF (i j k : State) : ℝ :=
  M i k * (1 + (q j / q i) * (q k / q j) / 2) -
  M j k * (1 + (q k / q j) / 2) - M i j * (1 - q j / q i)

def tripletG (i j k : State) : ℝ := M i k * (q j / q i) - M j k

/-- The original inclusive kappa branches and strict sign tests. -/
def originalTriplet (i j k : State) : Prop :=
  let κ := (M i k + M j k) / (M i j * (1 + q j / q i))
  (1 ≤ κ ∧ 0 < tripletF i j k) ∨ (κ ≤ 1 ∧ 0 < tripletG i j k)

theorem q_positive (i : State) : 0 < q i := by
  fin_cases i <;> norm_num [q]

theorem q_strictly_decreasing (i j : State) (h : i < j) : q j < q i := by
  fin_cases i <;> fin_cases j <;> norm_num [q, Fin.lt_def] at *

theorem positive_rates (i j : State) (h : i ≠ j) : 0 < M i j := by
  fin_cases i <;> fin_cases j <;> norm_num [M] at *

theorem conservation (j : State) : (∑ i, M i j) = 0 := by
  fin_cases j <;> norm_num [M, Fin.sum_univ_succ]

theorem detailed_balance (i j : State) : M i j * q j = M j i * q i := by
  fin_cases i <;> fin_cases j <;> norm_num [M, q]

theorem stationary (i : State) : (∑ j, M i j * q j) = 0 := by
  fin_cases i <;> norm_num [M, q, Fin.sum_univ_succ]

theorem triplet_signs (i j k : State) (hij : i < j) (hjk : j < k) :
    tripletF i j k < 0 ∧ tripletG i j k < 0 := by
  fin_cases i <;> fin_cases j <;> fin_cases k <;>
    norm_num [tripletF, tripletG, M, q, Fin.lt_def] at *

theorem all_original_triplets_fail (i j k : State) (hij : i < j) (hjk : j < k) :
    ¬ originalTriplet i j k := by
  obtain ⟨hf, hg⟩ := triplet_signs i j k hij hjk
  rintro (⟨_, h⟩ | ⟨_, h⟩)
  · exact (not_lt_of_ge hf.le) h
  · exact (not_lt_of_ge hg.le) h

theorem left_eigenvector (j : State) :
    (∑ i, M i j * w i) = slowRate * w j := by
  fin_cases j <;> norm_num [M, w, slowRate, Fin.sum_univ_succ]

theorem zero_bath_mean : (∑ i, q i * w i) = 0 := by
  norm_num [q, w, Fin.sum_univ_succ]

theorem mode_norm_identity : (∑ i, q i * w i ^ 2) = modeNorm := by
  norm_num [q, w, modeNorm, Fin.sum_univ_succ]

theorem right_eigenvector (i : State) :
    (∑ j, M i j * rightMode j) = slowRate * rightMode i := by
  fin_cases i <;> norm_num [M, q, w, rightMode, modeNorm, slowRate, Fin.sum_univ_succ]

theorem left_right_normalization : overlap rightMode = 1 := by
  norm_num [overlap, q, w, rightMode, modeNorm, Fin.sum_univ_succ]

theorem slowRate_location : -20 < slowRate ∧ slowRate < -19 := by
  norm_num [slowRate]

end
end Mpemba
