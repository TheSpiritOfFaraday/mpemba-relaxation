import Mpemba.Dynamics
import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Analysis.Calculus.Deriv.Prod
noncomputable section
namespace Mpemba
noncomputable def actLinear : (State → ℝ) →ₗ[ℝ] (State → ℝ) where
  toFun := act
  map_add' p r := by ext i; exact act_add p r i
  map_smul' c p := by ext i; exact act_mul c p i
noncomputable def actContinuous : (State → ℝ) →L[ℝ] (State → ℝ) :=
  actLinear.toContinuousLinearMap
theorem evolve_unique (p : State → ℝ) (f : ℝ → State → ℝ)
    (hf : ∀ t, HasDerivAt f (act (f t)) t) (h0 : f 0 = p) : f = evolve p := by
  apply ODE_solution_unique_univ (v := fun _ => act) (s := fun _ => Set.univ)
    (K := ‖actContinuous‖₊) (t₀ := 0)
  · intro t
    exact actContinuous.lipschitz.lipschitzOnWith
  · intro t
    exact ⟨hf t, Set.mem_univ _⟩
  · intro t
    refine ⟨hasDerivAt_pi.mpr (fun i => evolve_solves_master_equation p t i), Set.mem_univ _⟩
  · rw [h0]
    ext i
    exact (evolve_initial p i).symm

end Mpemba
