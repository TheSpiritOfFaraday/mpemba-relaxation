# Completion audit

All three requested formalization items are implemented and verified.
The final `./check.sh` run builds the complete library and passes 43 distinct
named-theorem axiom audits. Only `propext`, `Classical.choice`, and `Quot.sound`
occur. All 27 proof modules are imported by `Mpemba.lean`. A source scan finds
no `sorry`, `admit`, custom `axiom`, `unsafe`, or `native_decide` constructs.

## Requirement-by-requirement evidence

| Requested requirement | Authoritative Lean evidence |
| --- | --- |
| Literal Appendix B rational witness | `Data.lean`; compared with the current manuscript's Appendix B |
| Every transition rate positive | `positive_rates` |
| Column conservation | `conservation` |
| Detailed balance | `detailed_balance` |
| Every original triplet predicate fails | `all_original_triplets_fail`; inclusive kappa branches and strict tests |
| Exact slow eigenvector and eigenvalue | `left_eigenvector`, `right_eigenvector`, `left_right_normalization` |
| Simple slow mode | `slow_algebraic_multiplicity_one`, `slow_mode_is_slowest`; exact characteristic factorization |
| Actual Gibbs overlap reversal | `gibbs_weight_physical`, `strict_thermal_overlap_reversal`; integer-power enclosures in Lean |
| Overlap ordering implies eventual relaxation-distance ordering | `overlap_implies_eventual_distance_order` |
| Actual master-equation evolution | `evolve_solves_master_equation`, `evolve_unique`, `master_equation_distance_reversal` |
| Hotter preparation initially farther away | `companion_distance_reversal`; independently proved initial bounds in `InitialDistance.lean` |
| Reversal at a finite positive time | `reversal_at_a_finite_positive_time` |
| Quarter upper bound from original premise | `full_generator_quarter_bound`, `full_generator_amplitude_advantage_below_quarter` |
| Valid sharpness family | `Sharp.eventually_family_attainable`, including positive rates, full triplet failures, simple slow mode and finite ordered cooling temperatures |
| Limiting advantage 1/4 | `Sharp.advantage_limit`; proved hot/cold overlap limits -1/30 and -2/45 |
| Sharp supremum, not attained | `sharp_quarter_isLUB`, `sharp_quarter_not_attained`, `exists_advantage_arbitrarily_close` |

## Sharpness construction and scope

`SharpFamily.lean` uses a simpler rational family than the manuscript's
uniformly gapped construction, retaining its bath weights and left mode and
fixing the downward rate a12=1. Its exact eigenvalue is -3. The F124 and G234
constraints are equalities, and every original triplet predicate fails.
Continuity of thirteen strictly positive checks at epsilon zero establishes
a full interval of valid positive epsilon. `PositiveMode.lean` derives
simplicity and the variational top eigenvalue from the positive contrast mode
and connecting positive edges. No numerical eigenvalue or sampled sign is
accepted as an assumption.

`SharpThermal.lean` chooses epsilon=exp(-L^2), hot=L^(-3), cold=L^(-1) for
sufficiently large finite L. It proves all thermal limits using real analysis
in Lean. `SharpQuarter.lean` defines the set of attainable ordinary cooling
amplitude advantages and proves that its least upper bound is 1/4.

The upper-bound premise explicitly uses the standard variational definition
of a simple largest eigenvalue on the bath-orthogonal subspace. The bath
inverse temperature is normalized to one. Mode shape and overlap signs are
derived, not assumed. `FourStatePhysics.lean` separately verifies the general
generator's conservation, detailed balance, positive off-diagonal rates and
real-exponential Gibbs interpretation.

The physical reversal is proved for Pearson chi-square / weighted Euclidean
ordering under the Appendix B generator. The manuscript's general KL bridge,
uniform spectral-gap estimate, all-eight-strict perturbed sharpness family,
and time-advance corollaries are outside this requested formalization's
coverage. No formal claim of those stronger results is made.

## Reproduction

Run `./check.sh`. Lean 4.19.0 and mathlib commit
`c44e0c8ee63ca166450922a373c7409c5d26b00b` are pinned; the manifest pins all
transitive packages. Workspace-local elan is in `.tools/elan` and is selected
automatically by the script. README records the interpreted mathlib cache
download fallback for the macOS dyld issue.

No paper text was changed. Formal verification does not establish publication
priority, external peer review, or fidelity of the original conjecture's
transcription; those still require scholarly assessment. No memory updates
were requested or made.
