import Mpemba.CoolingOverlap

noncomputable section
open Matrix
namespace Mpemba.Quarter

def bathMean (q : FourWeights) (w : Fin 4 → ℝ) : ℝ :=
  w 0+q.t*w 1+q.u*w 2+q.v*w 3

def stateNorm (q : FourWeights) (w : Fin 4 → ℝ) : ℝ :=
  w ⬝ᵥ diagonal q.vector *ᵥ w

def stateRayleigh (q : FourWeights) (r : FourRates) (w : Fin 4 → ℝ) : ℝ :=
  w ⬝ᵥ diagonal q.vector *ᵥ ((generator q r)ᵀ*ᵥw)

/-- The usual variational definition of a simple slow left eigenmode on
    the bath-orthogonal subspace of the full four-state generator. -/
structure SimpleSlowMode (q : FourWeights) (r : FourRates) (lam : ℝ)
    (w : Fin 4 → ℝ) : Prop where
  nonzero : w ≠ 0
  mean_zero : bathMean q w = 0
  eigen : (generator q r)ᵀ*ᵥw = lam • w
  largest : ∀ z, bathMean q z = 0 → stateRayleigh q r z ≤ lam*stateNorm q z
  simple : ∀ z, bathMean q z = 0 → (generator q r)ᵀ*ᵥz = lam • z →
    ∃ k : ℝ, z = k • w

def modeCoordinates (q : FourWeights) (w : Fin 4 → ℝ) : Fin 3 → ℝ :=
  ![-q.t*(w 1+q.u/q.Q2*w 2+q.v/q.Q2*w 3),
    -q.u/q.Q2*w 2-q.u*q.v/(q.Q2*q.Q3)*w 3, -q.v/q.Q3*w 3]

theorem helmert_mean_zero (q : FourWeights) (hq : q.Valid) (x : Fin 3 → ℝ) :
    bathMean q (helmert q*ᵥx) = 0 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  simp [bathMean, helmert, Matrix.mulVec, dotProduct, Fin.sum_univ_three]
  field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv]
  unfold FourWeights.Q2 FourWeights.Q3
  ring

theorem helmert_injective (q : FourWeights) (hq : q.Valid) :
    Function.Injective (fun x : Fin 3 → ℝ => helmert q*ᵥx) := by
  intro x y hxy
  change helmert q*ᵥx = helmert q*ᵥy at hxy
  have hh : diagonal (eta q)*ᵥx = diagonal (eta q)*ᵥy := by
    rw [← helmert_gram q hq]
    simp only [← mulVec_mulVec]
    rw [hxy]
  ext i
  have hi := congrFun hh i
  simp only [mulVec_diagonal] at hi
  exact mul_left_cancel₀ (ne_of_gt (eta_positive q hq i)) hi

theorem modeCoordinates_reconstruct (q : FourWeights) (hq : q.Valid) (w : Fin 4 → ℝ)
    (hm : bathMean q w = 0) : helmert q*ᵥmodeCoordinates q w = w := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hw0 : w 0 = -q.t*w 1-q.u*w 2-q.v*w 3 := by
    unfold bathMean at hm
    linarith
  ext i
  fin_cases i <;> simp [helmert, modeCoordinates, Matrix.mulVec, dotProduct,
    Fin.sum_univ_three, hw0]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv,
    ne_of_gt hQ2, ne_of_gt hQ3]
  all_goals dsimp only [FourWeights.Q2, FourWeights.Q3]
  all_goals ring

theorem helmert_stateNorm (q : FourWeights) (hq : q.Valid) (x : Fin 3 → ℝ) :
    stateNorm q (helmert q*ᵥx) = weightedNorm q x := by
  have hid : x ⬝ᵥ diagonal (eta q)*ᵥx = stateNorm q (helmert q*ᵥx) := by
    rw [← helmert_gram q hq]
    simp only [← mulVec_mulVec]
    rw [dotProduct_mulVec, vecMul_transpose]
    rfl
  rw [← hid]
  simp only [weightedNorm, dotProduct, mulVec_diagonal]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem weightedDefect_hermitian (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (lam : ℝ) : (weightedDefect (eta q) (contrastMatrix q r) lam).IsHermitian := by
  apply Matrix.IsHermitian.ext
  intro i j
  simp only [star_trivial, weightedDefect_apply]
  by_cases hij : i = j
  · subst j
    rfl
  · rw [if_neg hij, if_neg (Ne.symm hij)]
    have hh := contrast_weighted_symmetric q hq r i j
    nlinarith only [hh]

theorem weightedDefect_quadratic (q : FourWeights) (r : FourRates) (lam : ℝ)
    (x : Fin 3 → ℝ) :
    x ⬝ᵥ weightedDefect (eta q) (contrastMatrix q r) lam *ᵥx =
      lam*weightedNorm q x-∑ i, eta q i*x i*(contrastMatrix q r*ᵥx) i := by
  simp only [dotProduct, weightedDefect_mulVec, weightedNorm,
    Finset.mul_sum, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem SimpleSlowMode.coordinates (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (lam : ℝ) (w : Fin 4 → ℝ) (h : SimpleSlowMode q r lam w) :
    SimpleTopMode (eta q) (contrastMatrix q r) lam (modeCoordinates q w) := by
  have hrec := modeCoordinates_reconstruct q hq w h.mean_zero
  have hinter (x : Fin 3 → ℝ) :
      (generator q r)ᵀ*ᵥ(helmert q*ᵥx) = helmert q*ᵥ(contrastMatrix q r*ᵥx) := by
    simp only [mulVec_mulVec]
    rw [helmert_intertwines q hq r]
  refine ⟨eta_positive q hq, ?_, ?_, ?_, ?_⟩
  · intro hz
    apply h.nonzero
    rw [← hrec, hz, mulVec_zero]
  · apply helmert_injective q hq
    change helmert q*ᵥ(contrastMatrix q r*ᵥmodeCoordinates q w) =
      helmert q*ᵥ(lam • modeCoordinates q w)
    rw [← hinter, hrec, h.eigen, mulVec_smul, hrec]
  · refine ⟨weightedDefect_hermitian q hq r lam, ?_⟩
    intro x
    simp only [star_trivial, weightedDefect_quadratic]
    have hh := h.largest (helmert q*ᵥx) (helmert_mean_zero q hq x)
    rw [helmert_stateNorm q hq x] at hh
    have hdir := contrast_dirichlet_identity q hq r x
    have hstate := generator_dirichlet_identity q hq r (helmert q*ᵥx)
    change -stateRayleigh q r (helmert q*ᵥx) = dirichletEnergy q r x at hstate
    linarith
  · intro y hy
    have hyeigen : (generator q r)ᵀ*ᵥ(helmert q*ᵥy) = lam • (helmert q*ᵥy) := by
      rw [hinter, hy, mulVec_smul]
    obtain ⟨k, hk⟩ := h.simple _ (helmert_mean_zero q hq y) hyeigen
    refine ⟨k, helmert_injective q hq ?_⟩
    change helmert q*ᵥy = helmert q*ᵥ(k • modeCoordinates q w)
    rw [mulVec_smul, hrec]
    exact hk

theorem SimpleTopMode.full_mode (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x) :
    SimpleSlowMode q r lam (helmert q*ᵥx) := by
  have hinter (y : Fin 3 → ℝ) :
      (generator q r)ᵀ*ᵥ(helmert q*ᵥy) = helmert q*ᵥ(contrastMatrix q r*ᵥy) := by
    simp only [mulVec_mulVec]
    rw [helmert_intertwines q hq r]
  refine ⟨?_, helmert_mean_zero q hq x, ?_, ?_, ?_⟩
  · intro hz
    apply h.nonzero
    apply helmert_injective q hq
    simpa only [mulVec_zero] using hz
  · rw [hinter, h.eigen, mulVec_smul]
  · intro z hz
    let y := modeCoordinates q z
    have hrec : helmert q*ᵥy = z := modeCoordinates_reconstruct q hq z hz
    have hquad := h.top.re_dotProduct_nonneg y
    simp only [star_trivial, RCLike.re_to_real, weightedDefect_quadratic] at hquad
    have hnorm := helmert_stateNorm q hq y
    have hstate := generator_dirichlet_identity q hq r (helmert q*ᵥy)
    change -stateRayleigh q r (helmert q*ᵥy) = dirichletEnergy q r y at hstate
    have hdir := contrast_dirichlet_identity q hq r y
    rw [hrec] at hstate hnorm
    rw [hnorm]
    linarith
  · intro z hz heigen
    let y := modeCoordinates q z
    have hrec : helmert q*ᵥy = z := modeCoordinates_reconstruct q hq z hz
    have hyeigen : contrastMatrix q r*ᵥy = lam • y := by
      apply helmert_injective q hq
      change helmert q*ᵥ(contrastMatrix q r*ᵥy) = helmert q*ᵥ(lam • y)
      rw [← hinter, hrec, heigen, mulVec_smul, hrec]
    obtain ⟨k, hk⟩ := h.simple y hyeigen
    refine ⟨k, ?_⟩
    rw [← hrec, hk, mulVec_smul]

/-- The upper bound stated for a full generator and its slow left
    eigenfunction, rather than supplied contrast coordinates. -/
theorem full_generator_quarter_bound (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Positive) (htrip : allOriginalTripletsFail q r) (lam : ℝ) (w : Fin 4 → ℝ)
    (h : SimpleSlowMode q r lam w) (hot cold : ℝ)
    (hh : 0 < hot) (hhc : hot ≤ cold) (hc : cold < 1) :
    |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) cold| /
      |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) hot| < 4/3 := by
  have hb := original_triplets_quarter_bound q hq r hr htrip lam (modeCoordinates q w)
    (h.coordinates q hq r lam w) hot cold hh hhc hc
  simpa only [modeOverlap, modeCoordinates_reconstruct q hq w h.mean_zero] using hb

theorem full_generator_amplitude_advantage_below_quarter
    (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Positive) (htrip : allOriginalTripletsFail q r) (lam : ℝ) (w : Fin 4 → ℝ)
    (h : SimpleSlowMode q r lam w) (hot cold : ℝ)
    (hh : 0 < hot) (hhc : hot ≤ cold) (hc : cold < 1) :
    1-|fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) hot| /
      |fullMean q.t q.u q.v (w 0) (w 1) (w 2) (w 3) cold| < 1/4 := by
  have htop := h.coordinates q hq r lam w
  have hcone := (allOriginalTripletsFail_iff_cone q hq r hr).mp htrip
  have hhot := simple_top_mode_overlap_nonzero q hq r hcone lam _ htop hot
    (lt_of_le_of_lt hhc hc)
  have hcold := simple_top_mode_overlap_nonzero q hq r hcone lam _ htop cold hc
  simp only [modeOverlap, modeCoordinates_reconstruct q hq w h.mean_zero] at hhot hcold
  exact Mpemba.amplitude_advantage_below_quarter _ _ (abs_pos.mpr hhot) (abs_pos.mpr hcold)
    (full_generator_quarter_bound q hq r hr htrip lam w h hot cold hh hhc hc)

end Mpemba.Quarter
