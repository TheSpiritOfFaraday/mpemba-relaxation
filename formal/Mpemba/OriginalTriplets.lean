import Mpemba.FourState

namespace Mpemba.Quarter

/-- The original inclusive kappa branches, before clearing any denominator. -/
def originalPredicate (q₁ q₂ q₃ A B C : ℝ) : Prop :=
  let t := q₂/q₁
  let u := q₃/q₂
  let F := A*(1+t*u/2)-B*(1+u/2)-C*(1-t)
  let G := A*t-B
  let κ := (A+B)/(C*(1+t))
  (1 ≤ κ ∧ 0 < F) ∨ (κ ≤ 1 ∧ 0 < G)

theorem original_failure_cleared (q₁ q₂ q₃ A B C : ℝ)
    (hq₁ : 0 < q₁) (hq₂ : 0 < q₂) (hq₃ : 0 < q₃) (h12 : q₂ < q₁)
    (hC : 0 < C) :
    ¬ originalPredicate q₁ q₂ q₃ A B C ↔
      A*q₂*(2*q₁+q₃)-B*q₁*(2*q₂+q₃)-2*C*q₂*(q₁-q₂) ≤ 0 ∧ A*q₂ ≤ B*q₁ := by
  have ht : 0 < q₂/q₁ := div_pos hq₂ hq₁
  have ht1 : q₂/q₁ < 1 := (div_lt_one hq₁).mpr h12
  have hu : 0 < q₃/q₂ := div_pos hq₃ hq₂
  have hF : 2*q₁*q₂*(A*(1+(q₂/q₁)*(q₃/q₂)/2)-B*(1+(q₃/q₂)/2)-C*(1-q₂/q₁)) =
      A*q₂*(2*q₁+q₃)-B*q₁*(2*q₂+q₃)-2*C*q₂*(q₁-q₂) := by
    field_simp [ne_of_gt hq₁, ne_of_gt hq₂]
    ring
  have hG : q₁*(A*(q₂/q₁)-B) = A*q₂-B*q₁ := by
    field_simp [ne_of_gt hq₁]
    ring
  have hiff := Mpemba.original_triplet_failure_iff A B C (q₂/q₁) (q₃/q₂) hC ht ht1 hu
  dsimp only at hiff
  dsimp only [originalPredicate]
  rw [hiff]
  constructor
  · rintro ⟨h₁, h₂⟩
    have hh₁ := mul_nonpos_of_nonneg_of_nonpos (by positivity : 0 ≤ 2*q₁*q₂) h₁
    have hh₂ := mul_nonpos_of_nonneg_of_nonpos hq₁.le h₂
    rw [hF] at hh₁
    rw [hG] at hh₂
    exact ⟨hh₁, by linarith⟩
  · rintro ⟨h₁, h₂⟩
    refine ⟨?_, ?_⟩
    · rw [← hF] at h₁
      exact nonpos_of_mul_nonpos_right h₁ (by positivity : 0 < 2*q₁*q₂)
    · have hh₂ : A*q₂-B*q₁ ≤ 0 := sub_nonpos.mpr h₂
      rw [← hG] at hh₂
      exact nonpos_of_mul_nonpos_right hh₂ hq₁

def FourRates.Positive (r : FourRates) : Prop :=
  0 < r.a ∧ 0 < r.b ∧ 0 < r.c ∧ 0 < r.d ∧ 0 < r.e ∧ 0 < r.f

theorem FourRates.Positive.nonnegative {r : FourRates} (h : r.Positive) : r.Nonnegative :=
  ⟨h.1.le, h.2.1.le, h.2.2.1.le, h.2.2.2.1.le, h.2.2.2.2.1.le, h.2.2.2.2.2.le⟩

def allOriginalTripletsFail (q : FourWeights) (r : FourRates) : Prop :=
  ¬ originalPredicate 1 q.t q.u r.b r.d r.a ∧
  ¬ originalPredicate 1 q.t q.v r.c r.e r.a ∧
  ¬ originalPredicate 1 q.u q.v r.c r.f r.b ∧
  ¬ originalPredicate q.t q.u q.v r.e r.f r.d

theorem allOriginalTripletsFail_iff_cone (q : FourWeights) (hq : q.Valid)
    (r : FourRates) (hr : r.Positive) : allOriginalTripletsFail q r ↔ TripletCone q r := by
  obtain ⟨ht, hu, hv, hQ2, hQ3, hQ4⟩ := hq.positive
  have ht1 := hq.2.2.2
  have hut := hq.2.2.1
  have hu1 := lt_trans hut ht1
  unfold allOriginalTripletsFail
  rw [original_failure_cleared 1 q.t q.u r.b r.d r.a (by norm_num) ht hu ht1 hr.1,
    original_failure_cleared 1 q.t q.v r.c r.e r.a (by norm_num) ht hv ht1 hr.1,
    original_failure_cleared 1 q.u q.v r.c r.f r.b (by norm_num) hu hv hu1 hr.2.1,
    original_failure_cleared q.t q.u q.v r.e r.f r.d ht hu hv hut hr.2.2.2.1]
  simp only [mul_one, one_mul]
  constructor
  · rintro ⟨⟨hF₁, hG₁⟩, ⟨hF₂, hG₂⟩, ⟨hF₃, hG₃⟩, ⟨hF₄, hG₄⟩⟩
    exact ⟨hG₁, hG₂, hG₃, hG₄, hF₁, hF₂, hF₃, hF₄⟩
  · intro h
    exact ⟨⟨h.f123, h.g123⟩, ⟨h.f124, h.g124⟩, ⟨h.f134, h.g134⟩, ⟨h.f234, h.g234⟩⟩

end Mpemba.Quarter
