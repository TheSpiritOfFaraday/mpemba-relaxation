import Mpemba.TopMode
import Mpemba.TripletLogic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity

noncomputable section
open Matrix
namespace Mpemba.Quarter

/-- Bath weights after the harmless common rescaling q₁=1. -/
structure FourWeights where
  t : ℝ
  u : ℝ
  v : ℝ

def FourWeights.Valid (q : FourWeights) : Prop := 0 < q.v ∧ q.v < q.u ∧ q.u < q.t ∧ q.t < 1
def FourWeights.vector (q : FourWeights) : Fin 4 → ℝ := ![1, q.t, q.u, q.v]
def FourWeights.Q2 (q : FourWeights) : ℝ := 1+q.t
def FourWeights.Q3 (q : FourWeights) : ℝ := 1+q.t+q.u
def FourWeights.Q4 (q : FourWeights) : ℝ := 1+q.t+q.u+q.v

theorem FourWeights.Valid.positive {q : FourWeights} (h : q.Valid) :
    0 < q.t ∧ 0 < q.u ∧ 0 < q.v ∧ 0 < q.Q2 ∧ 0 < q.Q3 ∧ 0 < q.Q4 := by
  unfold FourWeights.Valid at h
  unfold FourWeights.Q2 FourWeights.Q3 FourWeights.Q4
  rcases h with ⟨hv, hvu, hut, ht1⟩
  constructor
  · linarith
  constructor
  · linarith
  exact ⟨hv, by linarith, by linarith, by linarith⟩

structure FourRates where
  a : ℝ
  b : ℝ
  c : ℝ
  d : ℝ
  e : ℝ
  f : ℝ

def FourRates.Nonnegative (r : FourRates) : Prop :=
  0 ≤ r.a ∧ 0 ≤ r.b ∧ 0 ≤ r.c ∧ 0 ≤ r.d ∧ 0 ≤ r.e ∧ 0 ≤ r.f

/-- The eight F/G inequalities with positive denominators cleared. -/
structure TripletCone (q : FourWeights) (r : FourRates) : Prop where
  g123 : r.b*q.t ≤ r.d
  g124 : r.c*q.t ≤ r.e
  g134 : r.c*q.u ≤ r.f
  g234 : r.e*q.u ≤ r.f*q.t
  f123 : r.b*q.t*(2+q.u)-r.d*(2*q.t+q.u)-2*r.a*q.t*(1-q.t) ≤ 0
  f124 : r.c*q.t*(2+q.v)-r.e*(2*q.t+q.v)-2*r.a*q.t*(1-q.t) ≤ 0
  f134 : r.c*q.u*(2+q.v)-r.f*(2*q.u+q.v)-2*r.b*q.u*(1-q.u) ≤ 0
  f234 : r.e*q.u*(2*q.t+q.v)-r.f*q.t*(2*q.u+q.v)-2*r.d*q.u*(q.t-q.u) ≤ 0

def generator (q : FourWeights) (r : FourRates) : Matrix (Fin 4) (Fin 4) ℝ := !![
  -r.a*q.t-r.b*q.u-r.c*q.v, r.a, r.b, r.c;
  r.a*q.t, -r.a-r.d*q.u/q.t-r.e*q.v/q.t, r.d, r.e;
  r.b*q.u, r.d*q.u/q.t, -r.b-r.d-r.f*q.v/q.u, r.f;
  r.c*q.v, r.e*q.v/q.t, r.f*q.v/q.u, -r.c-r.e-r.f]

def helmert (q : FourWeights) : Matrix (Fin 4) (Fin 3) ℝ := !![
  1, 1, 1; -1/q.t, 1, 1; 0, -q.Q2/q.u, 1; 0, 0, -q.Q3/q.v]

def eta (q : FourWeights) : Fin 3 → ℝ :=
  ![q.Q2/q.t, q.Q2*q.Q3/q.u, q.Q3*q.Q4/q.v]

/-- Coordinate generator in the unnormalized weighted Helmert basis. -/
def contrastMatrix (q : FourWeights) (r : FourRates) : Matrix (Fin 3) (Fin 3) ℝ := !![
  -r.a*q.Q2-(r.b*q.t^2*q.u+r.d*q.u+r.c*q.t^2*q.v+r.e*q.v)/(q.t*q.Q2),
    (q.Q3*(r.d-r.b*q.t)+q.v*(r.e-r.c*q.t))/q.Q2,
    q.Q4*(r.e-r.c*q.t)/q.Q2;
  q.u*(q.Q3*(r.d-r.b*q.t)+q.v*(r.e-r.c*q.t))/(q.t*q.Q2*q.Q3),
    -(r.b+r.d)*q.Q3/q.Q2-
      q.v*(q.u^2*(r.c+r.e)+r.f*q.Q2^2)/(q.u*q.Q2*q.Q3),
    q.Q4*(r.f*q.Q2-q.u*(r.c+r.e))/(q.Q2*q.Q3);
  q.v*(r.e-r.c*q.t)/(q.t*q.Q3),
    q.v*(r.f*q.Q2-q.u*(r.c+r.e))/(q.u*q.Q3),
    -(r.c+r.e+r.f)*q.Q4/q.Q3]

theorem eta_positive (q : FourWeights) (hq : q.Valid) (i : Fin 3) : 0 < eta q i := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  fin_cases i <;> simp only [eta, cons_val_zero', cons_val_succ'] <;> positivity

set_option maxHeartbeats 0 in
theorem helmert_gram (q : FourWeights) (hq : q.Valid) :
    (helmert q)ᵀ * diagonal q.vector * helmert q = diagonal (eta q) := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_succ, helmert, FourWeights.vector, eta,
      Matrix.diagonal, FourWeights.Q2, FourWeights.Q3, FourWeights.Q4]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv] <;> ring

set_option maxHeartbeats 0 in
theorem helmert_intertwines (q : FourWeights) (hq : q.Valid) (r : FourRates) :
    (generator q r)ᵀ * helmert q = helmert q * contrastMatrix q r := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_succ, generator, helmert, contrastMatrix,
      FourWeights.Q2, FourWeights.Q3, FourWeights.Q4]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv,
    ne_of_gt hQ2, ne_of_gt hQ3]
  all_goals ring

set_option maxHeartbeats 0 in
theorem contrast_weighted_symmetric (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (i j : Fin 3) : eta q i * contrastMatrix q r i j = eta q j * contrastMatrix q r j i := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  fin_cases i <;> fin_cases j <;>
    simp [eta, contrastMatrix, FourWeights.Q2, FourWeights.Q3, FourWeights.Q4]
  all_goals field_simp [ne_of_gt ht, ne_of_gt hu, ne_of_gt hv,
    ne_of_gt hQ2, ne_of_gt hQ3]
  all_goals ring

theorem contrast_offdiagonal_nonnegative (q : FourWeights) (hq : q.Valid) (r : FourRates)
    (hcone : TripletCone q r) (i j : Fin 3) (hij : i ≠ j) :
    0 ≤ contrastMatrix q r i j := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have h₁ : 0 ≤ r.d-r.b*q.t := sub_nonneg.mpr hcone.g123
  have h₂ : 0 ≤ r.e-r.c*q.t := sub_nonneg.mpr hcone.g124
  have h₃ : 0 ≤ r.f*q.Q2-q.u*(r.c+r.e) := by
    unfold FourWeights.Q2
    nlinarith only [hcone.g134, hcone.g234]
  fin_cases i <;> fin_cases j <;> norm_num at hij <;>
    simp only [contrastMatrix, Matrix.of_apply, cons_val_zero', cons_val_succ'] <;> positivity

end Mpemba.Quarter
