import Mpemba.SlowModeBridge

noncomputable section
open Matrix
namespace Mpemba.Quarter

def threeKernelEnergy (A : Matrix (Fin 3) (Fin 3) ℝ) (x y : Fin 3 → ℝ) : ℝ :=
  (-A 0 1)*x 0*x 1*(y 0/x 0-y 1/x 1)^2+
    (-A 0 2)*x 0*x 2*(y 0/x 0-y 2/x 2)^2+
    (-A 1 2)*x 1*x 2*(y 1/x 1-y 2/x 2)^2

theorem three_kernel_energy_identity (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : A.IsHermitian) (x : Fin 3 → ℝ) (hx : ∀ i, x i ≠ 0) (hkernel : A*ᵥx = 0)
    (y : Fin 3 → ℝ) : y ⬝ᵥ A*ᵥy = threeKernelEnergy A x y := by
  have h₀ := congrFun hkernel 0
  have h₁ := congrFun hkernel 1
  have h₂ := congrFun hkernel 2
  simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three, Pi.zero_apply] at h₀ h₁ h₂
  have hs₁ : A 1 0 = A 0 1 := by simpa only [star_trivial] using hA.apply 0 1
  have hs₂ : A 2 0 = A 0 2 := by simpa only [star_trivial] using hA.apply 0 2
  have hs₃ : A 2 1 = A 1 2 := by simpa only [star_trivial] using hA.apply 1 2
  have hd₀ : A 0 0 = -(A 0 1*x 1+A 0 2*x 2)/x 0 := by
    apply (eq_div_iff (hx 0)).mpr
    linarith
  have hd₁ : A 1 1 = -(A 0 1*x 0+A 1 2*x 2)/x 1 := by
    apply (eq_div_iff (hx 1)).mpr
    rw [hs₁] at h₁
    linarith
  have hd₂ : A 2 2 = -(A 0 2*x 0+A 1 2*x 1)/x 2 := by
    apply (eq_div_iff (hx 2)).mpr
    rw [hs₂, hs₃] at h₂
    linarith
  simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_three, threeKernelEnergy]
  rw [hd₀, hd₁, hd₂, hs₁, hs₂, hs₃]
  field_simp [hx 0, hx 1, hx 2]
  ring

theorem positive_kernel_psd_simple (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : A.IsHermitian) (hoff : ∀ i j, i ≠ j → A i j ≤ 0)
    (h₀₁ : A 0 1 < 0) (h₀₂ : A 0 2 < 0)
    (x : Fin 3 → ℝ) (hx : ∀ i, 0 < x i) (hkernel : A*ᵥx = 0) :
    A.PosSemidef ∧ ∀ y, A*ᵥy = 0 → ∃ k : ℝ, y = k • x := by
  have hid := three_kernel_energy_identity A hA x (fun i => ne_of_gt (hx i)) hkernel
  have hcoef₁ : 0 < (-A 0 1)*x 0*x 1 :=
    mul_pos (mul_pos (neg_pos.mpr h₀₁) (hx 0)) (hx 1)
  have hcoef₂ : 0 < (-A 0 2)*x 0*x 2 :=
    mul_pos (mul_pos (neg_pos.mpr h₀₂) (hx 0)) (hx 2)
  have hcoef₃ : 0 ≤ (-A 1 2)*x 1*x 2 := by
    exact mul_nonneg (mul_nonneg (neg_nonneg.mpr (hoff 1 2 (by decide))) (hx 1).le) (hx 2).le
  have hnonneg (y : Fin 3 → ℝ) : 0 ≤ threeKernelEnergy A x y := by
    unfold threeKernelEnergy
    positivity
  refine ⟨⟨hA, ?_⟩, ?_⟩
  · intro y
    simp only [star_trivial, hid]
    exact hnonneg y
  · intro y hy
    have hzero : threeKernelEnergy A x y = 0 := by
      rw [← hid y, hy]
      simp
    have hn₁ := mul_nonneg hcoef₁.le (sq_nonneg (y 0/x 0-y 1/x 1))
    have hn₂ := mul_nonneg hcoef₂.le (sq_nonneg (y 0/x 0-y 2/x 2))
    have hn₃ := mul_nonneg hcoef₃ (sq_nonneg (y 1/x 1-y 2/x 2))
    unfold threeKernelEnergy at hzero
    have he₁ : (-A 0 1)*x 0*x 1*(y 0/x 0-y 1/x 1)^2 = 0 := by linarith
    have he₂ : (-A 0 2)*x 0*x 2*(y 0/x 0-y 2/x 2)^2 = 0 := by linarith
    have hd₁ : y 0/x 0 = y 1/x 1 := by
      have hh := (mul_eq_zero.mp he₁).resolve_left (ne_of_gt hcoef₁)
      exact sub_eq_zero.mp (sq_eq_zero_iff.mp hh)
    have hd₂ : y 0/x 0 = y 2/x 2 := by
      have hh := (mul_eq_zero.mp he₂).resolve_left (ne_of_gt hcoef₂)
      exact sub_eq_zero.mp (sq_eq_zero_iff.mp hh)
    refine ⟨y 0/x 0, ?_⟩
    ext i
    fin_cases i
    · change y 0 = (y 0/x 0)*x 0
      exact (div_mul_cancel₀ _ (ne_of_gt (hx 0))).symm
    · change y 1 = (y 0/x 0)*x 1
      rw [hd₁, div_mul_cancel₀ _ (ne_of_gt (hx 1))]
    · change y 2 = (y 0/x 0)*x 2
      rw [hd₂, div_mul_cancel₀ _ (ne_of_gt (hx 2))]

/-- A positive contrast eigenvector, with two positive edges connecting the
    three coordinates, certifies both the top eigenvalue and simplicity. -/
theorem simpleTop_of_positive_mode (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (lam : ℝ) (x : Fin 3 → ℝ) (hx : ∀ i, 0 < x i)
    (heigen : contrastMatrix q r*ᵥx = lam • x)
    (hoff : ∀ i j, i ≠ j → 0 ≤ contrastMatrix q r i j)
    (h₀₁ : 0 < contrastMatrix q r 0 1) (h₀₂ : 0 < contrastMatrix q r 0 2) :
    SimpleTopMode (eta q) (contrastMatrix q r) lam x := by
  let A := weightedDefect (eta q) (contrastMatrix q r) lam
  have hkernel : A*ᵥx = 0 :=
    (weightedDefect_kernel_iff (eta q) (eta_positive q hq) _ lam x).mpr heigen
  have hAoff (i j : Fin 3) (hij : i ≠ j) : A i j ≤ 0 := by
    dsimp only [A]
    rw [weightedDefect_apply, if_neg hij, zero_sub]
    exact mul_nonpos_of_nonneg_of_nonpos (eta_positive q hq i).le (neg_nonpos.mpr (hoff i j hij))
  have hA₀₁ : A 0 1 < 0 := by
    dsimp only [A]
    rw [weightedDefect_apply, if_neg (by decide : (0:Fin 3) ≠ 1), zero_sub]
    exact mul_neg_of_pos_of_neg (eta_positive q hq 0) (neg_neg_of_pos h₀₁)
  have hA₀₂ : A 0 2 < 0 := by
    dsimp only [A]
    rw [weightedDefect_apply, if_neg (by decide : (0:Fin 3) ≠ 2), zero_sub]
    exact mul_neg_of_pos_of_neg (eta_positive q hq 0) (neg_neg_of_pos h₀₂)
  obtain ⟨htop, hsimple⟩ := positive_kernel_psd_simple A
    (weightedDefect_hermitian q hq r lam) hAoff hA₀₁ hA₀₂ x hx hkernel
  refine ⟨eta_positive q hq, ?_, heigen, htop, ?_⟩
  · intro hz
    have hh := congrFun hz 0
    exact (ne_of_gt (hx 0)) hh
  · intro y hy
    exact hsimple y ((weightedDefect_kernel_iff (eta q) (eta_positive q hq) _ lam y).mpr hy)

end Mpemba.Quarter
