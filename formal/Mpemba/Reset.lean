import Mpemba.FourState

noncomputable section
open Matrix
namespace Mpemba.Quarter

def FourRates.reset (r : FourRates) (q : FourWeights) (k : ℝ) : FourRates :=
  ⟨r.a-k, r.b-k, r.c-k, r.d-k*q.t, r.e-k*q.t, r.f-k*q.u⟩

def FourRates.resetMinimum (r : FourRates) : ℝ := min r.a (min r.b r.c)

theorem reset_preserves_cone (q : FourWeights) (r : FourRates) (k : ℝ)
    (h : TripletCone q r) : TripletCone q (r.reset q k) := by
  constructor <;> dsimp only [FourRates.reset]
  · nlinarith only [h.g123]
  · nlinarith only [h.g124]
  · nlinarith only [h.g134]
  · nlinarith only [h.g234]
  · nlinarith only [h.f123]
  · nlinarith only [h.f124]
  · nlinarith only [h.f134]
  · nlinarith only [h.f234]

theorem reset_minimum_bounds (r : FourRates) :
    r.resetMinimum ≤ r.a ∧ r.resetMinimum ≤ r.b ∧ r.resetMinimum ≤ r.c :=
  ⟨min_le_left _ _, le_trans (min_le_right _ _) (min_le_left _ _),
    le_trans (min_le_right _ _) (min_le_right _ _)⟩

theorem reset_minimum_nonnegative (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (_hr : r.Nonnegative) (hcone : TripletCone q r) :
    (r.reset q r.resetMinimum).Nonnegative := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨ha, hb, hc⟩ := reset_minimum_bounds r
  unfold FourRates.Nonnegative FourRates.reset
  refine ⟨sub_nonneg.mpr ha, sub_nonneg.mpr hb, sub_nonneg.mpr hc, ?_, ?_, ?_⟩
  · exact sub_nonneg.mpr (le_trans (mul_le_mul_of_nonneg_right hb ht.le) hcone.g123)
  · exact sub_nonneg.mpr (le_trans (mul_le_mul_of_nonneg_right hc ht.le) hcone.g124)
  · exact sub_nonneg.mpr (le_trans (mul_le_mul_of_nonneg_right hc hu.le) hcone.g134)

theorem reset_has_zero_first_row (q : FourWeights) (r : FourRates) :
    (r.reset q r.resetMinimum).a = 0 ∨
    (r.reset q r.resetMinimum).b = 0 ∨ (r.reset q r.resetMinimum).c = 0 := by
  unfold FourRates.reset FourRates.resetMinimum
  by_cases ha : r.a ≤ min r.b r.c
  · simp [min_eq_left ha]
  · rw [min_eq_right (le_of_not_ge ha)]
    by_cases hb : r.b ≤ r.c
    · simp [min_eq_left hb]
    · simp [min_eq_right (le_of_not_ge hb)]

set_option maxHeartbeats 0 in
theorem contrast_reset (q : FourWeights) (hq : q.Valid) (r : FourRates) (k : ℝ) :
    contrastMatrix q (r.reset q k) = contrastMatrix q r + diagonal (fun _ => k*q.Q4) := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [contrastMatrix, FourRates.reset, Matrix.diagonal,
      FourWeights.Q2, FourWeights.Q3, FourWeights.Q4]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv,
    ne_of_gt hQ2, ne_of_gt hQ3]
  all_goals ring_nf
  all_goals simp

theorem weightedDefect_shift {n : Type*} [Fintype n] [DecidableEq n]
    (η : n → ℝ) (B : Matrix n n ℝ) (lam k : ℝ) :
    weightedDefect η (B+diagonal (fun _ => k)) (lam+k) = weightedDefect η B lam := by
  ext i j
  simp only [weightedDefect_apply, add_apply, diagonal_apply]
  split_ifs <;> ring

theorem SimpleTopMode.shift {n : Type*} [Fintype n] [DecidableEq n]
    {η : n → ℝ} {B : Matrix n n ℝ} {lam : ℝ} {x : n → ℝ}
    (h : SimpleTopMode η B lam x) (k : ℝ) :
    SimpleTopMode η (B+diagonal (fun _ => k)) (lam+k) x := by
  refine ⟨h.weights_pos, h.nonzero, ?_, ?_, ?_⟩
  · rw [add_mulVec, diagonal_const_mulVec, h.eigen, add_smul]
  · simpa only [weightedDefect_shift] using h.top
  · intro y hy
    apply h.simple
    ext i
    have hi := congrFun hy i
    simp only [add_mulVec, diagonal_const_mulVec, Pi.add_apply, Pi.smul_apply,
      smul_eq_mul] at hi ⊢
    linarith

theorem simpleTop_reset (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (lam : ℝ) (x : Fin 3 → ℝ) (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (k : ℝ) : SimpleTopMode (eta q) (contrastMatrix q (r.reset q k)) (lam+k*q.Q4) x := by
  rw [contrast_reset q hq r k]
  exact h.shift _

def weightedNorm (q : FourWeights) (x : Fin 3 → ℝ) : ℝ := ∑ i, eta q i * (x i)^2

def stateDirichletEnergy (q : FourWeights) (r : FourRates) (w : Fin 4 → ℝ) : ℝ :=
  r.a*q.t*(w 0-w 1)^2 + r.b*q.u*(w 0-w 2)^2 + r.c*q.v*(w 0-w 3)^2 +
    r.d*q.u*(w 1-w 2)^2 + r.e*q.v*(w 1-w 3)^2 + r.f*q.v*(w 2-w 3)^2

def dirichletEnergy (q : FourWeights) (r : FourRates) (x : Fin 3 → ℝ) : ℝ :=
  stateDirichletEnergy q r (helmert q *ᵥ x)

set_option maxHeartbeats 0 in
theorem generator_dirichlet_identity (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (w : Fin 4 → ℝ) :
    -(w ⬝ᵥ diagonal q.vector *ᵥ ((generator q r)ᵀ *ᵥ w)) = stateDirichletEnergy q r w := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  simp [stateDirichletEnergy, Matrix.mulVec, dotProduct, Fin.sum_univ_four, generator,
    FourWeights.vector, Matrix.diagonal]
  field_simp [ne_of_gt ht, ne_of_gt hu]
  ring

theorem contrast_dirichlet_identity (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (x : Fin 3 → ℝ) :
    -(∑ i, eta q i * x i * (contrastMatrix q r *ᵥ x) i) = dirichletEnergy q r x := by
  have hmat : diagonal (eta q)*contrastMatrix q r =
      (helmert q)ᵀ*diagonal q.vector*(generator q r)ᵀ*helmert q := by
    rw [← helmert_gram q hq]
    calc
      _ = ((helmert q)ᵀ*diagonal q.vector)*(helmert q*contrastMatrix q r) := by
        simp only [Matrix.mul_assoc]
      _ = ((helmert q)ᵀ*diagonal q.vector)*((generator q r)ᵀ*helmert q) := by
        rw [helmert_intertwines q hq r]
      _ = _ := by simp only [Matrix.mul_assoc]
  have hsum : (∑ i, eta q i * x i * (contrastMatrix q r *ᵥ x) i) =
      x ⬝ᵥ (diagonal (eta q)*contrastMatrix q r) *ᵥ x := by
    rw [← mulVec_mulVec]
    simp only [dotProduct, mulVec_diagonal]
    apply Finset.sum_congr rfl
    intro i _
    ring
  have hform : x ⬝ᵥ (diagonal (eta q)*contrastMatrix q r) *ᵥ x =
      (helmert q*ᵥx) ⬝ᵥ diagonal q.vector *ᵥ ((generator q r)ᵀ *ᵥ (helmert q*ᵥx)) := by
    rw [hmat]
    simp only [← mulVec_mulVec]
    rw [dotProduct_mulVec, vecMul_transpose]
  rw [hsum, hform]
  exact generator_dirichlet_identity q hq r _

theorem dirichletEnergy_nonnegative (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (x : Fin 3 → ℝ) : 0 ≤ dirichletEnergy q r x := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨ha, hb, hc, hd, he, hf⟩ := hr
  unfold dirichletEnergy stateDirichletEnergy
  positivity

theorem weightedNorm_positive (q : FourWeights) (hq : q.Valid) (x : Fin 3 → ℝ)
    (hx : x ≠ 0) : 0 < weightedNorm q x := by
  have hex : ∃ i, x i ≠ 0 := by
    by_contra hn
    apply hx
    ext i
    by_contra hi
    exact hn ⟨i, hi⟩
  obtain ⟨i, hi⟩ := hex
  have hterm := Finset.single_le_sum (f := fun j => eta q j * (x j)^2)
    (fun j _ => mul_nonneg (eta_positive q hq j).le (sq_nonneg _)) (Finset.mem_univ i)
  exact lt_of_lt_of_le (mul_pos (eta_positive q hq i) (sq_pos_of_ne_zero hi)) hterm

theorem contrast_eigenvalue_nonpositive (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (lam : ℝ) (x : Fin 3 → ℝ) (hx : x ≠ 0)
    (heigen : contrastMatrix q r *ᵥ x = lam • x) : lam ≤ 0 := by
  have heq : (∑ i, eta q i * x i * (contrastMatrix q r *ᵥ x) i) = lam*weightedNorm q x := by
    rw [heigen]
    simp only [weightedNorm, Finset.mul_sum, Pi.smul_apply, smul_eq_mul]
    apply Finset.sum_congr rfl
    intro i _
    ring
  have hdir := contrast_dirichlet_identity q hq r x
  rw [heq] at hdir
  have hnonneg := dirichletEnergy_nonnegative q hq r hr x
  have hn := weightedNorm_positive q hq x hx
  by_contra hl
  have hp := mul_pos (lt_of_not_ge hl) hn
  linarith

end Mpemba.Quarter
