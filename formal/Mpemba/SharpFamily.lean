import Mpemba.SlowModeBridge
import Mpemba.PositiveMode
import Mathlib.Tactic.FunProp

noncomputable section
open Matrix Filter Topology
namespace Mpemba.Quarter.Sharp

/-- A simpler rational family with the manuscript's limiting Gibbs geometry.
    The fixed downward rate a=1 leaves two exact boundary triplet equalities. -/
def weights (ε : ℝ) : FourWeights := ⟨1-ε/3, 3/5, ε^2⟩
def mode (ε : ℝ) : Fin 4 → ℝ := ![1, -4/5, -1/3-4*ε/9+5*ε^3/3, -ε]

def P (ε : ℝ) : ℝ := 15*ε^4-5*ε^3+69*ε^2-20*ε+68
def R (ε : ℝ) : ℝ := 15*ε^3-4*ε-12
def S (ε : ℝ) : ℝ := 75*ε^3-20*ε+21

def rates (ε : ℝ) : FourRates := ⟨1,
  (150*ε^6+640*ε^5+190*ε^4+459*ε^3+138*ε^2-252*ε-1224)/(R ε*P ε),
  (30*ε^3+125*ε^2-42*ε+240)/(3*P ε),
  5*(3-ε)*(225*ε^5-185*ε^4+425*ε^3-113*ε^2-60*ε+204)/(S ε*P ε),
  5*(3-ε)*(9*ε^2-2*ε+16)/(3*P ε), 3*(9*ε^2-2*ε+16)/P ε⟩

theorem zero_rates : rates 0 = ⟨1,3/2,20/17,15/7,20/17,12/17⟩ := by
  norm_num [rates, P, R, S]

theorem weights_valid (ε : ℝ) (hε : 0 < ε) (hsmall : ε < 1/30) : (weights ε).Valid := by
  dsimp only [weights, FourWeights.Valid]
  have hsq : ε^2 < 1/900 := by nlinarith
  exact ⟨sq_pos_of_pos hε, by linarith, by linarith, by linarith⟩

theorem mode_mean_zero (ε : ℝ) : bathMean (weights ε) (mode ε) = 0 := by
  simp [bathMean, weights, mode]
  ring

theorem mode_nonzero (ε : ℝ) : mode ε ≠ 0 := by
  intro h
  have hh := congrFun h 0
  norm_num [mode] at hh

set_option maxHeartbeats 0 in
theorem family_eigen (ε : ℝ) (hP : P ε ≠ 0) (hR : R ε ≠ 0) (hS : S ε ≠ 0)
    (ht : 1-ε/3 ≠ 0) :
    (generator (weights ε) (rates ε))ᵀ*ᵥmode ε = (-3:ℝ) • mode ε := by
  have ht' : -ε/3+1 ≠ 0 := by intro h; apply ht; linarith
  have h3 : 3-ε ≠ 0 := by intro h; apply ht; linarith
  have h3' : -ε+3 ≠ 0 := by intro h; apply ht; linarith
  ext i
  fin_cases i <;> simp [generator, weights, rates, mode, Matrix.mulVec, dotProduct,
    Fin.sum_univ_four]
  all_goals field_simp [hP, hR, hS, ht, ht', h3, h3']
  all_goals dsimp only [P, R, S]
  all_goals ring

theorem boundary_F124 (ε : ℝ) (hP : P ε ≠ 0) :
    (rates ε).c*(weights ε).t*(2+(weights ε).v)-
      (rates ε).e*(2*(weights ε).t+(weights ε).v)-
      2*(rates ε).a*(weights ε).t*(1-(weights ε).t) = 0 := by
  simp [rates, weights]
  field_simp [hP]
  unfold P
  ring

theorem boundary_G234 (ε : ℝ) :
    (rates ε).e*(weights ε).u = (rates ε).f*(weights ε).t := by
  simp [rates, weights]
  ring

theorem G124_factor (ε : ℝ) (hP : P ε ≠ 0) :
    (rates ε).e-(rates ε).c*(weights ε).t =
      ε*(2*(ε-3)*(15*ε^2-5*ε-6)/(9*P ε)) := by
  simp [rates, weights]
  field_simp [hP]
  ring

def F123 (ε : ℝ) : ℝ :=
  (rates ε).b*(weights ε).t*(2+(weights ε).u)-
    (rates ε).d*(2*(weights ε).t+(weights ε).u)-
    2*(rates ε).a*(weights ε).t*(1-(weights ε).t)
def F134 (ε : ℝ) : ℝ :=
  (rates ε).c*(weights ε).u*(2+(weights ε).v)-
    (rates ε).f*(2*(weights ε).u+(weights ε).v)-
    2*(rates ε).b*(weights ε).u*(1-(weights ε).u)
def F234 (ε : ℝ) : ℝ :=
  (rates ε).e*(weights ε).u*(2*(weights ε).t+(weights ε).v)-
    (rates ε).f*(weights ε).t*(2*(weights ε).u+(weights ε).v)-
    2*(rates ε).d*(weights ε).u*((weights ε).t-(weights ε).u)

theorem zero_triplet_values : F123 0 = -117/70 ∧ F134 0 = -66/425 ∧ F234 0 = -276/595 := by
  norm_num [F123, F134, F234, zero_rates, weights]

theorem family_continuous_at_zero :
    ContinuousAt (fun ε => (rates ε).b) 0 ∧ ContinuousAt (fun ε => (rates ε).c) 0 ∧
    ContinuousAt (fun ε => (rates ε).d) 0 ∧ ContinuousAt (fun ε => (rates ε).e) 0 ∧
    ContinuousAt (fun ε => (rates ε).f) 0 := by
  unfold rates P R S
  constructor
  · fun_prop (disch := norm_num)
  constructor
  · fun_prop (disch := norm_num)
  constructor
  · fun_prop (disch := norm_num)
  constructor <;> fun_prop (disch := norm_num)

/-- Each entry is strictly positive at zero. This yields an open interval
    certificate by continuity, rather than a finite parameter mesh. -/
def checks (ε : ℝ) : Fin 13 → ℝ := ![P ε, -R ε, S ε,
  (rates ε).b, (rates ε).c, (rates ε).d, (rates ε).e, (rates ε).f,
  (rates ε).d-(rates ε).b*(weights ε).t, -F123 ε, -F134 ε, -F234 ε,
  2*(ε-3)*(15*ε^2-5*ε-6)/(9*P ε)]

theorem checks_positive_zero (i : Fin 13) : 0 < checks 0 i := by
  fin_cases i <;> norm_num [checks, rates, weights, F123, F134, F234, P, R, S]

theorem checks_continuous_at_zero (i : Fin 13) : ContinuousAt (fun ε => checks ε i) 0 := by
  fin_cases i <;> simp [checks, rates, weights, F123, F134, F234, P, R, S] <;>
    fun_prop (disch := norm_num)

theorem eventually_checks_positive : ∀ᶠ ε in 𝓝 (0:ℝ), ∀ i, 0 < checks ε i := by
  apply Filter.eventually_all.mpr
  intro i
  exact (tendsto_order.mp (checks_continuous_at_zero i)).1 0 (checks_positive_zero i)

theorem rates_positive_of_checks (ε : ℝ) (h : ∀ i, 0 < checks ε i) : (rates ε).Positive := by
  exact ⟨by norm_num [rates], h 3, h 4, h 5, h 6, h 7⟩

theorem cone_of_checks (ε : ℝ) (hε : 0 < ε) (hsmall : ε < 1/30)
    (h : ∀ i, 0 < checks ε i) : TripletCone (weights ε) (rates ε) := by
  have hq := weights_valid ε hε hsmall
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hP : P ε ≠ 0 := ne_of_gt (h 0)
  have h124 : 0 < (rates ε).e-(rates ε).c*(weights ε).t := by
    rw [G124_factor ε hP]
    exact mul_pos hε (h 12)
  have hG234 := boundary_G234 ε
  have h134 : (rates ε).c*(weights ε).u < (rates ε).f := by
    have hg := mul_pos h124 hu
    apply (mul_lt_mul_right ht).mp
    nlinarith only [hg, hG234]
  constructor
  · have hh := h 8
    change 0 < (rates ε).d-(rates ε).b*(weights ε).t at hh
    linarith
  · linarith
  · exact h134.le
  · exact hG234.le
  · exact (neg_pos.mp (h 9)).le
  · exact (boundary_F124 ε hP).le
  · exact (neg_pos.mp (h 10)).le
  · exact (neg_pos.mp (h 11)).le

theorem contrast_first_edges_positive (ε : ℝ) (hε : 0 < ε) (hsmall : ε < 1/30)
    (h : ∀ i, 0 < checks ε i) :
    0 < contrastMatrix (weights ε) (rates ε) 0 1 ∧
    0 < contrastMatrix (weights ε) (rates ε) 0 2 := by
  have hq := weights_valid ε hε hsmall
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hg₁ : 0 < (rates ε).d-(rates ε).b*(weights ε).t := h 8
  have hg₂ : 0 < (rates ε).e-(rates ε).c*(weights ε).t := by
    rw [G124_factor ε (ne_of_gt (h 0))]
    exact mul_pos hε (h 12)
  constructor <;> simp only [contrastMatrix, Matrix.of_apply]
  · change 0 < ((weights ε).Q3*((rates ε).d-(rates ε).b*(weights ε).t)+
      (weights ε).v*((rates ε).e-(rates ε).c*(weights ε).t))/(weights ε).Q2
    positivity
  · change 0 < (weights ε).Q4*((rates ε).e-(rates ε).c*(weights ε).t)/(weights ε).Q2
    positivity

theorem coordinates_positive (ε : ℝ) (hε : 0 < ε) (hsmall : ε < 1/30) :
    ∀ i, 0 < modeCoordinates (weights ε) (mode ε) i := by
  have hq := weights_valid ε hε hsmall
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have hε1 : ε < 1 := by linarith
  have hsq : ε^2 < 1 := by nlinarith
  have hcube : ε^3 < ε := by
    have hh := mul_lt_mul_of_pos_right hsq hε
    nlinarith only [hh]
  have hr : mode ε 2 < 0 := by
    change -1/3-4*ε/9+5*ε^3/3 < 0
    nlinarith only [hcube, hsmall]
  have hlast : mode ε 3 < 0 := by change -ε < 0; linarith
  have hrec := modeCoordinates_reconstruct (weights ε) hq (mode ε) (mode_mean_zero ε)
  have hdiff := (mode_shape_identities (weights ε) hq (modeCoordinates (weights ε) (mode ε))).1
  rw [hrec] at hdiff
  have hx₀ : 0 < modeCoordinates (weights ε) (mode ε) 0 := by
    have hprod : 0 < (weights ε).Q2/(weights ε).t*modeCoordinates (weights ε) (mode ε) 0 := by
      rw [← hdiff]
      norm_num [mode]
    exact (mul_pos_iff_of_pos_left (div_pos hQ2 ht)).mp hprod
  intro i
  fin_cases i
  · exact hx₀
  · change 0 < -(weights ε).u/(weights ε).Q2*mode ε 2-
      (weights ε).u*(weights ε).v/((weights ε).Q2*(weights ε).Q3)*mode ε 3
    have hfirst := mul_pos_of_neg_of_neg (div_neg_of_neg_of_pos (neg_neg_of_pos hu) hQ2) hr
    have hsecond : (weights ε).u*(weights ε).v/((weights ε).Q2*(weights ε).Q3)*mode ε 3 ≤ 0 :=
      mul_nonpos_of_nonneg_of_nonpos (by positivity) hlast.le
    linarith
  · change 0 < -(weights ε).v/(weights ε).Q3*mode ε 3
    exact mul_pos_of_neg_of_neg (div_neg_of_neg_of_pos (neg_neg_of_pos hv) hQ3) hlast

theorem family_simple_slow (ε : ℝ) (hε : 0 < ε) (hsmall : ε < 1/30)
    (h : ∀ i, 0 < checks ε i) : SimpleSlowMode (weights ε) (rates ε) (-3) (mode ε) := by
  have hq := weights_valid ε hε hsmall
  have ht : 1-ε/3 ≠ 0 := ne_of_gt hq.positive.1
  have heigen := family_eigen ε (ne_of_gt (h 0)) (ne_of_lt (neg_pos.mp (h 1)))
    (ne_of_gt (h 2)) ht
  have hrec := modeCoordinates_reconstruct (weights ε) hq (mode ε) (mode_mean_zero ε)
  have hcoord : contrastMatrix (weights ε) (rates ε)*ᵥmodeCoordinates (weights ε) (mode ε) =
      (-3:ℝ) • modeCoordinates (weights ε) (mode ε) := by
    apply helmert_injective (weights ε) hq
    change helmert (weights ε)*ᵥ(contrastMatrix (weights ε) (rates ε)*ᵥ_) =
      helmert (weights ε)*ᵥ((-3:ℝ) • modeCoordinates (weights ε) (mode ε))
    rw [mulVec_mulVec, ← helmert_intertwines (weights ε) hq (rates ε),
      ← mulVec_mulVec, hrec, heigen, mulVec_smul, hrec]
  obtain ⟨hedge₁, hedge₂⟩ := contrast_first_edges_positive ε hε hsmall h
  have htop := simpleTop_of_positive_mode (weights ε) hq (rates ε) (-3)
    (modeCoordinates (weights ε) (mode ε)) (coordinates_positive ε hε hsmall) hcoord
    (contrast_offdiagonal_nonnegative (weights ε) hq (rates ε) (cone_of_checks ε hε hsmall h))
    hedge₁ hedge₂
  have hfull := htop.full_mode (weights ε) hq (rates ε) (-3) _
  rwa [hrec] at hfull

end Mpemba.Quarter.Sharp
