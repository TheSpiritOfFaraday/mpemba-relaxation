import Mpemba.Data
import Mathlib.LinearAlgebra.Matrix.Charpoly.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
noncomputable section
namespace Mpemba
set_option maxHeartbeats 0
theorem characteristic_factorization (x : ℝ) : (Matrix.diagonal (fun _ : State => x) - M).det =
    x * (250000*x+4966671) *
    (2530028648200000000*x^2+1002356769715761942000*x+50981800562345565497767) /
    632507162050000000000000 := by
  have heq : Matrix.diagonal (fun _ : State => x) - M = !![
    x+49029/4000, -1, -28, -249/10;
    -199/200, x+10689318350801/117041850000, -330154675693/1470375000, -4983/200;
    -56/5, -660309351386/7315115625, x+36129268114499561/143029257750000, -1527327066/151990625;
    -249/4000, -4983/79600, -763663533/12159250000, x+72789920403/1215925000] := by
    ext i j
    fin_cases i <;> fin_cases j <;> norm_num [Matrix.diagonal, M]
  rw [heq, Matrix.det_succ_row_zero]
  simp [Matrix.det_fin_three, Fin.sum_univ_succ, Matrix.submatrix_apply,
    Matrix.diagonal, M, Fin.succAbove]
  ring

/-- Quadratic factor for the two fast eigenvalues. -/
def fastQuadratic (x : ℝ) : ℝ :=
  2530028648200000000*x^2 + 1002356769715761942000*x + 50981800562345565497767

theorem fast_quadratic_positive_above (x : ℝ) (hx : -59 ≤ x) :
    0 < fastQuadratic x := by
  have ht : 0 ≤ x + 59 := by linarith
  have hs : 0 ≤ (x + 59)^2 := sq_nonneg _
  unfold fastQuadratic
  nlinarith

/-- Every real characteristic root other than stationarity and the rational
slow eigenvalue is strictly faster, with a gap exceeding 39. -/
theorem other_characteristic_roots_faster (x : ℝ)
    (hroot : (Matrix.diagonal (fun _ : State => x) - M).det = 0)
    (hzero : x ≠ 0) (hslow : x ≠ slowRate) : x < -59 := by
  rw [characteristic_factorization] at hroot
  have hx : 250000*x+4966671 ≠ 0 := by
    intro h
    apply hslow
    unfold slowRate
    linarith
  have hr : fastQuadratic x = 0 := by
    have h := (div_eq_zero_iff).mp hroot
    rcases h with h | h
    · rcases mul_eq_zero.mp h with h | h
      · rcases mul_eq_zero.mp h with h | h
        · exact False.elim (hzero h)
        · exact False.elim (hx h)
      · exact h
    · norm_num at h
  by_contra h
  have hp := fast_quadratic_positive_above x (le_of_not_gt h)
  linarith

/-- The full left eigenspace is exactly the line through the printed mode. -/
theorem slow_eigenspace_one_dimensional (f : State → ℝ)
    (hf : ∀ j, (∑ i, M i j * f i) = slowRate * f j) :
    ∀ i, f i = f 0 * w i := by
  have h0 := hf 0
  have h1 := hf 1
  have h2 := hf 2
  norm_num [M, slowRate, Fin.sum_univ_succ, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.vecHead, Matrix.vecTail] at h0 h1 h2
  simp only [show (Fin.succ (2 : Fin 3) : State) = 3 from rfl] at h0 h1 h2
  intro i
  fin_cases i
  · simp [w]
  · change f 1 = f 0 * (-759/1000)
    linarith
  · change f 2 = f 0 * (-48957/80000)
    linarith
  · change f 3 = f 0 * (-1/250)
    linarith

end Mpemba
