import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Matrix.Notation
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

noncomputable section
open Matrix
namespace Mpemba.Quarter

/-- A PSD matrix with nonpositive off-diagonal entries has a kernel stable
    under coordinatewise absolute value. This is the Perron step needed here. -/
theorem psd_kernel_abs {n : Type*} [Fintype n] [DecidableEq n]
    (A : Matrix n n ℝ) (hA : A.PosSemidef)
    (hoff : ∀ i j, i ≠ j → A i j ≤ 0) (x : n → ℝ) (hx : A *ᵥ x = 0) :
    A *ᵥ (fun i => |x i|) = 0 := by
  have hterm (i j : n) : |x i| * (A i j * |x j|) ≤ x i * (A i j * x j) := by
    by_cases hij : i = j
    · subst j
      have hs := sq_abs (x i)
      nlinarith only [congrArg (fun z : ℝ => A i i * z) hs]
    · have h := mul_le_mul_of_nonpos_right (le_abs_self (x i * x j)) (hoff i j hij)
      rw [abs_mul] at h
      nlinarith only [h]
  have hquad : (fun i => |x i|) ⬝ᵥ A *ᵥ (fun i => |x i|) ≤ x ⬝ᵥ A *ᵥ x := by
    simp only [dotProduct, mulVec, Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i _
    apply Finset.sum_le_sum
    intro j _
    exact hterm i j
  have hzero : x ⬝ᵥ A *ᵥ x = 0 := by rw [hx]; simp
  rw [hzero] at hquad
  have hnonneg : 0 ≤ (fun i => |x i|) ⬝ᵥ A *ᵥ (fun i => |x i|) := by
    simpa only [RCLike.re_to_real] using hA.re_dotProduct_nonneg (fun i => |x i|)
  have hz := le_antisymm hquad hnonneg
  apply (hA.dotProduct_mulVec_zero_iff _).mp
  simpa only [star_trivial] using hz

def weightedDefect {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (B : Matrix n n ℝ) (lam : ℝ) : Matrix n n ℝ :=
  diagonal η * (diagonal (fun _ => lam) - B)

theorem weightedDefect_mulVec {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (B : Matrix n n ℝ) (lam : ℝ) (x : n → ℝ) (i : n) :
    (weightedDefect η B lam *ᵥ x) i = η i * (lam*x i - (B *ᵥ x) i) := by
  rw [weightedDefect, ← mulVec_mulVec, mulVec_diagonal, sub_mulVec]
  simp [mulVec_diagonal]

theorem weightedDefect_kernel_iff {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (hη : ∀ i, 0 < η i) (B : Matrix n n ℝ) (lam : ℝ) (x : n → ℝ) :
    weightedDefect η B lam *ᵥ x = 0 ↔ B *ᵥ x = lam • x := by
  constructor
  · intro h
    ext i
    have hi := congrFun h i
    rw [weightedDefect_mulVec] at hi
    have hz := (mul_eq_zero.mp hi).resolve_left (ne_of_gt (hη i))
    change (B *ᵥ x) i = lam*x i
    linarith
  · intro h
    ext i
    rw [weightedDefect_mulVec, h]
    simp

/-- A simple top eigenmode stated by its exact variational certificate.
    `top` is the PSD form of the largest-eigenvalue condition, not a sign or
    component-shape assumption about the mode. -/
structure SimpleTopMode {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (B : Matrix n n ℝ) (lam : ℝ) (x : n → ℝ) : Prop where
  weights_pos : ∀ i, 0 < η i
  nonzero : x ≠ 0
  eigen : B *ᵥ x = lam • x
  top : (weightedDefect η B lam).PosSemidef
  simple : ∀ y, B *ᵥ y = lam • y → ∃ c : ℝ, y = c • x

theorem SimpleTopMode.abs_eigen {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) (hoff : ∀ i j, i ≠ j → 0 ≤ B i j) :
    B *ᵥ (fun i => |x i|) = lam • (fun i => |x i|) := by
  apply (weightedDefect_kernel_iff η h.weights_pos B lam _).mp
  apply psd_kernel_abs _ h.top
  · intro i j hij
    simp only [weightedDefect, diagonal_mul, sub_apply, diagonal_apply_ne _ hij,
      zero_sub]
    exact mul_nonpos_of_nonneg_of_nonpos (h.weights_pos i).le (neg_nonpos.mpr (hoff i j hij))
  · exact (weightedDefect_kernel_iff η h.weights_pos B lam x).mpr h.eigen

theorem SimpleTopMode.abs_nonzero {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) : (fun i => |x i|) ≠ 0 := by
  intro hz
  apply h.nonzero
  ext i
  have hi := congrFun hz i
  exact abs_eq_zero.mp hi

theorem SimpleTopMode.abs_is_nonzero_multiple {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) (hoff : ∀ i j, i ≠ j → 0 ≤ B i j) :
    ∃ c : ℝ, c ≠ 0 ∧ (fun i => |x i|) = c • x := by
  obtain ⟨c, hc⟩ := h.simple _ (h.abs_eigen hoff)
  refine ⟨c, ?_, hc⟩
  intro hz
  apply h.abs_nonzero
  simpa only [hz, zero_smul] using hc

theorem SimpleTopMode.abs_mode {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) (hoff : ∀ i j, i ≠ j → 0 ≤ B i j) :
    SimpleTopMode η B lam (fun i => |x i|) := by
  obtain ⟨k, hk, hscale⟩ := h.abs_is_nonzero_multiple hoff
  refine ⟨h.weights_pos, h.abs_nonzero, h.abs_eigen hoff, h.top, ?_⟩
  intro y hy
  obtain ⟨c, hc⟩ := h.simple y hy
  refine ⟨c/k, ?_⟩
  rw [hc, hscale]
  simp only [smul_smul, div_mul_cancel₀ _ hk]

theorem psd_diagonal_nonneg {n : Type*} [Fintype n] [DecidableEq n]
    (A : Matrix n n ℝ) (hA : A.PosSemidef) (i : n) : 0 ≤ A i i := by
  have h := hA.re_dotProduct_nonneg (Pi.single i 1)
  simpa [star_trivial, single_dotProduct, mulVec_single] using h

theorem psd_zero_diagonal_kernel {n : Type*} [Fintype n] [DecidableEq n]
    (A : Matrix n n ℝ) (hA : A.PosSemidef) (i : n) (hi : A i i = 0) :
    A *ᵥ Pi.single i 1 = 0 := by
  apply (hA.dotProduct_mulVec_zero_iff _).mp
  simpa [star_trivial, single_dotProduct, mulVec_single] using hi

theorem weightedDefect_apply {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (B : Matrix n n ℝ) (lam : ℝ) (i j : n) :
    weightedDefect η B lam i j = η i * ((if i=j then lam else 0)-B i j) := by
  simp [weightedDefect, diagonal_mul, diagonal_apply]

/-- Simplicity makes the diagonal Rayleigh bound strict whenever the mode
    has a nonzero component outside that diagonal coordinate. -/
theorem SimpleTopMode.diagonal_gap {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) (j k : n) (hkj : k ≠ j) (hxk : x k ≠ 0) :
    0 < lam-B j j := by
  have hd := psd_diagonal_nonneg _ h.top j
  rw [weightedDefect_apply, if_pos rfl] at hd
  have hnonneg : 0 ≤ lam-B j j := by
    by_contra hn
    have hh := mul_neg_of_pos_of_neg (h.weights_pos j) (lt_of_not_ge hn)
    linarith
  apply lt_of_le_of_ne hnonneg
  intro heq
  have hdiag : weightedDefect η B lam j j = 0 := by
    rw [weightedDefect_apply, if_pos rfl]
    rw [← heq, mul_zero]
  have hkernel := psd_zero_diagonal_kernel _ h.top j hdiag
  have heigen := (weightedDefect_kernel_iff η h.weights_pos B lam _).mp hkernel
  obtain ⟨c, hc⟩ := h.simple _ heigen
  have hck := congrFun hc k
  simp [Pi.single_apply, hkj, Ne.symm hkj] at hck
  have hc0 : c = 0 := hck.resolve_right hxk
  have hcj := congrFun hc j
  simp [hc0] at hcj

end Mpemba.Quarter
