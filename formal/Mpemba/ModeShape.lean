import Mpemba.Reset
import Mpemba.SchurSpectral
import Mpemba.QuarterBoundary

noncomputable section
open Matrix
namespace Mpemba.Quarter

def shapeSlope (q : FourWeights) : ℝ := (1-q.t)*q.u/(2*q.t*q.Q3)

def shapeC (q : FourWeights) (r : FourRates) : ℝ :=
  let B := contrastMatrix q r
  B 1 0+shapeSlope q*(B 1 1-B 0 0)-(shapeSlope q)^2*B 0 1

def shapeD (q : FourWeights) (r : FourRates) : ℝ :=
  let B := contrastMatrix q r
  (B 1 2-shapeSlope q*B 0 2)*(B 2 0+shapeSlope q*B 2 1)

theorem shapeSlope_positive (q : FourWeights) (hq : q.Valid) : 0 < shapeSlope q := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hdiff : 0 < 1-q.t := sub_pos.mpr hq.2.2.2
  unfold shapeSlope
  positivity

set_option maxHeartbeats 0 in
theorem shapeD_first_factor (q : FourWeights) (hq : q.Valid) (r : FourRates) :
    contrastMatrix q r 1 2-shapeSlope q*contrastMatrix q r 0 2 =
      q.Q4/(2*q.t*q.Q3)*(2*r.f*q.t-r.c*q.t*q.u-r.e*q.u) := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  simp [contrastMatrix, shapeSlope, FourWeights.Q2, FourWeights.Q3, FourWeights.Q4]
  field_simp [ne_of_gt ht, ne_of_gt hQ2, ne_of_gt hQ3]
  ring

theorem shapeD_nonnegative (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hcone : TripletCone q r) : 0 ≤ shapeD q r := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hbracket : 0 ≤ 2*r.f*q.t-r.c*q.t*q.u-r.e*q.u := by
    have hg := mul_le_mul_of_nonneg_right hcone.g134 ht.le
    nlinarith only [hg, hcone.g234]
  dsimp only [shapeD]
  rw [shapeD_first_factor q hq r]
  exact mul_nonneg (mul_nonneg (by positivity) hbracket)
    (add_nonneg (contrast_offdiagonal_nonnegative q hq r hcone 2 0 (by decide))
      (mul_nonneg (shapeSlope_positive q hq).le
        (contrast_offdiagonal_nonnegative q hq r hcone 2 1 (by decide))))

def shapeK (q : FourWeights) (r : FourRates) : ℝ :=
  q.u*(r.e*(2+q.u)-r.c*q.t*(2*q.t+q.u))-2*r.f*q.t*(1-q.t)

def shapeK1 (q : FourWeights) (r : FourRates) : ℝ :=
  2*r.f*q.t-r.c*q.t*q.u-r.e*q.u

def shapeK2 (q : FourWeights) (r : FourRates) : ℝ :=
  2*q.Q3*(r.e-r.c*q.t)+(1-q.t)*(r.f*q.Q2-q.u*(r.c+r.e))

set_option maxHeartbeats 0 in
theorem shapeC_formula (q : FourWeights) (hq : q.Valid) (r : FourRates) :
    shapeC q r = q.Q2/(4*q.t^2*q.Q3^2)*
      (q.Q3*boundaryN 1 q.t q.u r.a r.b r.d+q.v*shapeK q r) := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  simp [shapeC, contrastMatrix, shapeSlope]
  field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hQ2, ne_of_gt hQ3]
  dsimp only [shapeK, boundaryN, FourWeights.Q2, FourWeights.Q3]
  ring

theorem shapeD_formula (q : FourWeights) (hq : q.Valid) (r : FourRates) :
    shapeD q r = q.Q4*q.v/(4*q.t^2*q.Q3^3)*shapeK1 q r*shapeK2 q r := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  dsimp only [shapeD]
  rw [shapeD_first_factor q hq r]
  simp [shapeSlope, contrastMatrix, shapeK1, shapeK2]
  field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hQ3]
  ring

theorem shapeK_identity (q : FourWeights) (r : FourRates) :
    (r.c+r.e+r.f)*q.Q2*shapeK q r+shapeK1 q r*shapeK2 q r =
      q.Q2*q.Q3*(2*r.c*r.e*q.u*(1-q.t)+r.f*boundaryJ 1 q.t q.u r.c r.e) := by
  unfold shapeK shapeK1 shapeK2 boundaryJ FourWeights.Q2 FourWeights.Q3
  ring

theorem shape_boundary_identity (q : FourWeights) (hq : q.Valid) (r : FourRates) :
    -contrastMatrix q r 2 2*shapeC q r+shapeD q r =
      q.Q2*q.Q4/(4*q.t^2*q.Q3^2)*
        boundaryT 1 q.t q.u q.v r.a r.b r.c r.d r.e r.f := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hb : -contrastMatrix q r 2 2 = (r.c+r.e+r.f)*q.Q4/q.Q3 := by
    simp [contrastMatrix]
    ring
  rw [hb, shapeC_formula q hq r, shapeD_formula q hq r]
  calc
    _ = q.Q4/(4*q.t^2*q.Q3^3)*
        ((r.c+r.e+r.f)*q.Q2*q.Q3*boundaryN 1 q.t q.u r.a r.b r.d+
          q.v*((r.c+r.e+r.f)*q.Q2*shapeK q r+shapeK1 q r*shapeK2 q r)) := by
      field_simp [ne_of_gt ht, ne_of_gt hQ3]
      ring
    _ = _ := by
      rw [shapeK_identity]
      unfold boundaryT
      field_simp [ne_of_gt ht, ne_of_gt hQ3]
      ring

theorem shape_boundary_nonnegative (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r)
    (hface : r.a = 0 ∨ r.b = 0 ∨ r.c = 0) :
    0 ≤ -contrastMatrix q r 2 2*shapeC q r+shapeD q r := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨ha, hb, hc, hd, he, hf⟩ := hr
  have hN : 0 ≤ boundaryN 1 q.t q.u r.a r.b r.d := by
    have h := mul_nonpos_of_nonpos_of_nonneg hcone.f123 hu.le
    unfold boundaryN
    nlinarith only [h]
  have hF : r.c*q.t*(2*1+q.v) ≤ r.e*1*(2*q.t+q.v)+2*r.a*q.t*(1-q.t) := by
    nlinarith only [hcone.f124]
  have hT := boundaryT_nonnegative (by norm_num : (0:ℝ)<1) ht hu hv
    hq.2.2.2.le hq.2.1.le ha hb hc hd he hf hN hF hface
  rw [shape_boundary_identity q hq r]
  exact mul_nonneg (by positivity) hT

theorem boundary_mode_ratio (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r)
    (hface : r.a = 0 ∨ r.b = 0 ∨ r.c = 0) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (hx : ∀ i, 0 ≤ x i) : shapeSlope q*x 0 ≤ x 1 := by
  by_cases hx0 : x 0 = 0
  · simpa only [hx0, mul_zero] using hx 1
  have hxpos := lt_of_le_of_ne (hx 0) (Ne.symm hx0)
  have hδ := h.diagonal_gap 2 0 (by decide) hx0
  have hlam := contrast_eigenvalue_nonpositive q hq r hr lam x h.nonzero h.eigen
  have hsign := schur_sign_from_boundary lam (contrastMatrix q r 2 2) (shapeC q r)
    (shapeD q r) hlam hδ.le (shapeD_nonnegative q hq r hcone)
    (shape_boundary_nonnegative q hq r hr hcone hface)
  apply schur_ratio_comparison (eta q) (contrastMatrix q r) lam x h
    (contrast_offdiagonal_nonnegative q hq r hcone) hxpos (shapeSlope q)
    (shapeSlope_positive q hq).le
  have hid := schur_ratio_numerator
    (contrastMatrix q r 0 0) (contrastMatrix q r 0 1) (contrastMatrix q r 0 2)
    (contrastMatrix q r 1 0) (contrastMatrix q r 1 1) (contrastMatrix q r 1 2)
    (contrastMatrix q r 2 0) (contrastMatrix q r 2 1)
    (lam-contrastMatrix q r 2 2) (shapeSlope q) (ne_of_gt hδ)
  dsimp only [shapeC, shapeD] at hsign
  rw [← hid] at hsign
  have hh := nonneg_of_mul_nonneg_right hsign hδ
  simpa only [schurEntry, Fin.castSucc_zero, Fin.castSucc_one] using hh

/-- The arithmetic restriction is derived from the generator inequalities,
    including every reducible reset boundary. -/
theorem nonnegative_mode_ratio (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (hx : ∀ i, 0 ≤ x i) : shapeSlope q*x 0 ≤ x 1 :=
  boundary_mode_ratio q hq (r.reset q r.resetMinimum)
    (reset_minimum_nonnegative q hq r hr hcone)
    (reset_preserves_cone q r r.resetMinimum hcone) (reset_has_zero_first_row q r)
    (lam+r.resetMinimum*q.Q4) x (simpleTop_reset q hq r lam x h r.resetMinimum) hx

theorem mode_shape_identities (q : FourWeights) (hq : q.Valid) (x : Fin 3 → ℝ) :
    (helmert q*ᵥx) 0-(helmert q*ᵥx) 1 = q.Q2/q.t*x 0 ∧
    (helmert q*ᵥx) 0+(helmert q*ᵥx) 1-2*(helmert q*ᵥx) 2 =
      2*q.Q3/q.u*(x 1-shapeSlope q*x 0) ∧
    (helmert q*ᵥx) 3 = -q.Q3/q.v*x 2 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  simp [helmert, Matrix.mulVec, dotProduct, Fin.sum_univ_three, shapeSlope,
    FourWeights.Q2, FourWeights.Q3]
  constructor
  · field_simp [ne_of_gt ht]; ring
  · field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hQ3]; ring

theorem nonnegative_mode_shape (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hr : r.Nonnegative) (hcone : TripletCone q r) (lam : ℝ) (x : Fin 3 → ℝ)
    (h : SimpleTopMode (eta q) (contrastMatrix q r) lam x)
    (hx : ∀ i, 0 ≤ x i) :
    (helmert q*ᵥx) 1 ≤ (helmert q*ᵥx) 0 ∧
    2*(helmert q*ᵥx) 2 ≤ (helmert q*ᵥx) 0+(helmert q*ᵥx) 1 ∧
    (helmert q*ᵥx) 3 ≤ 0 := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  obtain ⟨h₁, h₂, h₃⟩ := mode_shape_identities q hq x
  have hg₁ : 0 ≤ q.Q2/q.t*x 0 := mul_nonneg (by positivity) (hx 0)
  have hg₂ : 0 ≤ 2*q.Q3/q.u*(x 1-shapeSlope q*x 0) :=
    mul_nonneg (by positivity) (sub_nonneg.mpr (nonnegative_mode_ratio q hq r hr hcone lam x h hx))
  have hg₃ : -q.Q3/q.v*x 2 ≤ 0 :=
    mul_nonpos_of_nonpos_of_nonneg
      (div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr hQ3.le) hv.le) (hx 2)
  exact ⟨by linarith, by linarith, by linarith⟩

end Mpemba.Quarter
