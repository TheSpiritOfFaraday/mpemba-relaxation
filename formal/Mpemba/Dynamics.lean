import Mpemba.Thermal
import Mpemba.Spectral
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Tactic.FieldSimp

noncomputable section
namespace Mpemba

set_option maxHeartbeats 0

def fastLinear : ℝ := 1002356769715761942000 / 2530028648200000000
def fastConstant : ℝ := 50981800562345565497767 / 2530028648200000000
def fastDiscriminant : ℝ := fastLinear ^ 2 - 4 * fastConstant
def fastRate₁ : ℝ := (-fastLinear - Real.sqrt fastDiscriminant) / 2
def fastRate₂ : ℝ := (-fastLinear + Real.sqrt fastDiscriminant) / 2

theorem fast_discriminant_positive : 0 < fastDiscriminant := by
  norm_num [fastDiscriminant, fastLinear, fastConstant]

theorem fast_sqrt_bounds : 276 < Real.sqrt fastDiscriminant ∧
    Real.sqrt fastDiscriminant < 277 := by
  have hs := Real.sq_sqrt fast_discriminant_positive.le
  have hn := Real.sqrt_nonneg fastDiscriminant
  have hl : (276 : ℝ)^2 < fastDiscriminant := by
    norm_num [fastDiscriminant, fastLinear, fastConstant]
  have hu : fastDiscriminant < (277 : ℝ)^2 := by
    norm_num [fastDiscriminant, fastLinear, fastConstant]
  constructor <;> nlinarith

theorem fast_rates_distinct : fastRate₁ ≠ fastRate₂ := by
  have hs := fast_sqrt_bounds.1
  unfold fastRate₁ fastRate₂
  linarith

theorem fast_rates_sum : fastRate₁ + fastRate₂ = -fastLinear := by
  unfold fastRate₁ fastRate₂
  ring

theorem fast_rates_product : fastRate₁ * fastRate₂ = fastConstant := by
  have hs := Real.sq_sqrt fast_discriminant_positive.le
  change Real.sqrt fastDiscriminant ^ 2 = fastLinear ^ 2 - 4 * fastConstant at hs
  unfold fastRate₁ fastRate₂
  nlinarith

theorem fast_rates_below : fastRate₁ < -59 ∧ fastRate₂ < -59 := by
  obtain ⟨hl, hu⟩ := fast_sqrt_bounds
  unfold fastRate₁ fastRate₂
  norm_num [fastLinear] at *
  constructor <;> linarith

def act (p : State → ℝ) (i : State) : ℝ := ∑ j, M i j * p j
def total (p : State → ℝ) : ℝ := ∑ i, p i
def bath (i : State) : ℝ := q i / (∑ j, q j)
def fastPart (p : State → ℝ) (i : State) : ℝ :=
  p i - total p * bath i - overlap p * rightMode i

/-- Exact annihilating quadratic on the complementary two-dimensional space. -/
theorem fast_part_quadratic (p : State → ℝ) (i : State) :
    act (act (fastPart p)) i + fastLinear * act (fastPart p) i +
      fastConstant * fastPart p i = 0 := by
  fin_cases i <;>
    norm_num [act, fastPart, total, bath, overlap, rightMode, modeNorm,
      M, q, w, fastLinear, fastConstant, Fin.sum_univ_succ] <;>
    simp only [show (Fin.succ (2 : Fin 3) : State) = 3 from rfl,
      show (⟨2, by decide⟩ : State) = 2 from rfl,
      show (⟨3, by decide⟩ : State) = 3 from rfl] <;> ring

def fastComponent₁ (p : State → ℝ) (i : State) : ℝ :=
  (act (fastPart p) i - fastRate₂ * fastPart p i) / (fastRate₁ - fastRate₂)

def fastComponent₂ (p : State → ℝ) (i : State) : ℝ :=
  (fastRate₁ * fastPart p i - act (fastPart p) i) / (fastRate₁ - fastRate₂)

theorem fast_components_sum (p : State → ℝ) (i : State) :
    fastComponent₁ p i + fastComponent₂ p i = fastPart p i := by
  unfold fastComponent₁ fastComponent₂
  field_simp [sub_ne_zero.mpr fast_rates_distinct]
  ring

theorem act_add (p r : State → ℝ) (i : State) :
    act (fun j => p j + r j) i = act p i + act r i := by
  simp [act, mul_add, Finset.sum_add_distrib]

theorem act_sub (p r : State → ℝ) (i : State) :
    act (fun j => p j - r j) i = act p i - act r i := by
  simp [act, mul_sub, Finset.sum_sub_distrib]

theorem act_mul (c : ℝ) (p : State → ℝ) (i : State) :
    act (fun j => c * p j) i = c * act p i := by
  simp only [act, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

theorem act_div (p : State → ℝ) (c : ℝ) (i : State) :
    act (fun j => p j / c) i = act p i / c := by
  simp [act, mul_div_assoc, Finset.sum_div]

theorem fast_component₁_eigen (p : State → ℝ) (i : State) :
    act (fastComponent₁ p) i = fastRate₁ * fastComponent₁ p i := by
  unfold fastComponent₁
  rw [act_div, act_sub, act_mul]
  have hq := fast_part_quadratic p i
  have he : act (act (fastPart p)) i =
      (fastRate₁ + fastRate₂) * act (fastPart p) i -
      (fastRate₁ * fastRate₂) * fastPart p i := by
    rw [fast_rates_sum, fast_rates_product]
    linarith
  rw [he]
  ring

theorem fast_component₂_eigen (p : State → ℝ) (i : State) :
    act (fastComponent₂ p) i = fastRate₂ * fastComponent₂ p i := by
  unfold fastComponent₂
  rw [act_div, act_sub, act_mul]
  have hq := fast_part_quadratic p i
  have he : act (act (fastPart p)) i =
      (fastRate₁ + fastRate₂) * act (fastPart p) i -
      (fastRate₁ * fastRate₂) * fastPart p i := by
    rw [fast_rates_sum, fast_rates_product]
    linarith
  rw [he]
  ring

/-- Explicit solution of the master equation; the following theorems establish
    the initial condition and generator equation. -/
def evolve (p : State → ℝ) (t : ℝ) (i : State) : ℝ :=
  total p * bath i + Real.exp (slowRate*t) * (overlap p * rightMode i) +
    Real.exp (fastRate₁*t) * fastComponent₁ p i +
    Real.exp (fastRate₂*t) * fastComponent₂ p i

theorem evolve_initial (p : State → ℝ) (i : State) : evolve p 0 i = p i := by
  simp only [evolve, mul_zero, Real.exp_zero, one_mul]
  have h := fast_components_sum p i
  unfold fastPart at h
  linarith

theorem bath_stationary (i : State) : act bath i = 0 := by
  unfold bath
  rw [act_div]
  change (∑ j, M i j * q j) / _ = 0
  rw [stationary, zero_div]

theorem evolve_generator (p : State → ℝ) (t : ℝ) (i : State) :
    act (evolve p t) i =
      slowRate * Real.exp (slowRate*t) * (overlap p * rightMode i) +
      fastRate₁ * Real.exp (fastRate₁*t) * fastComponent₁ p i +
      fastRate₂ * Real.exp (fastRate₂*t) * fastComponent₂ p i := by
  unfold evolve
  rw [act_add, act_add, act_add, act_mul, act_mul, act_mul, act_mul, act_mul]
  rw [bath_stationary, fast_component₁_eigen, fast_component₂_eigen]
  have hr : act rightMode i = slowRate * rightMode i := right_eigenvector i
  rw [hr]
  ring

theorem evolve_solves_master_equation (p : State → ℝ) (t : ℝ) (i : State) :
    HasDerivAt (fun s => evolve p s i) (act (evolve p t) i) t := by
  rw [evolve_generator]
  have hs := (((hasDerivAt_id t).const_mul slowRate).exp).mul_const
    (overlap p * rightMode i)
  have h₁ := (((hasDerivAt_id t).const_mul fastRate₁).exp).mul_const (fastComponent₁ p i)
  have h₂ := (((hasDerivAt_id t).const_mul fastRate₂).exp).mul_const (fastComponent₂ p i)
  convert (((hasDerivAt_const t (total p * bath i)).add hs).add h₁).add h₂ using 1
  simp only [evolve, mul_one, zero_add, id_eq]
  ring

end Mpemba
