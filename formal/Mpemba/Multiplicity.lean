import Mpemba.Spectral
import Mathlib.LinearAlgebra.Matrix.Charpoly.Coeff
import Mathlib.Algebra.Polynomial.FieldDivision

noncomputable section
open Polynomial
namespace Mpemba

def certifiedCharpoly : ℝ[X] :=
  X * (C 250000 * X + C 4966671) *
    (C 2530028648200000000 * X^2 + C 1002356769715761942000 * X +
      C 50981800562345565497767) * C (1/632507162050000000000000)

theorem charpoly_eval_determinant (x : ℝ) :
    M.charpoly.eval x = (Matrix.diagonal (fun _ : State => x) - M).det := by
  rw [Matrix.charpoly, Matrix.eval_det, Matrix.matPolyEquiv_charmatrix]
  simp only [Polynomial.eval_sub, Polynomial.eval_X, Polynomial.eval_C, Matrix.scalar_apply]

theorem charpoly_exact : M.charpoly = certifiedCharpoly := by
  apply Polynomial.funext
  intro x
  rw [charpoly_eval_determinant, characteristic_factorization]
  simp [certifiedCharpoly]
  ring

theorem slow_charpoly_root : M.charpoly.IsRoot slowRate := by
  norm_num [Polynomial.IsRoot, charpoly_exact, certifiedCharpoly, slowRate]

theorem slow_charpoly_derivative_nonzero : M.charpoly.derivative.eval slowRate ≠ 0 := by
  norm_num [charpoly_exact, certifiedCharpoly, slowRate, Polynomial.derivative_mul,
    Polynomial.derivative_add, Polynomial.derivative_pow]

/-- Algebraic multiplicity, stronger than just uniqueness of the eigenvector line. -/
theorem slow_algebraic_multiplicity_one : M.charpoly.rootMultiplicity slowRate = 1 := by
  have hn : M.charpoly ≠ 0 := (Matrix.charpoly_monic M).ne_zero
  have hp := (Polynomial.rootMultiplicity_pos hn).mpr slow_charpoly_root
  have hnot : ¬ 1 < M.charpoly.rootMultiplicity slowRate := by
    intro h
    have hd := ((Polynomial.one_lt_rootMultiplicity_iff_isRoot hn).mp h).2
    exact slow_charpoly_derivative_nonzero hd
  omega

end Mpemba
