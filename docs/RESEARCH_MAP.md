# Research map

[Project home](../README.md) · [Documentation hub](README.md) · [Reproduction guide](../REPRODUCING.md)

Start with the [focused paper](../output/pdf/focused_paper.pdf) for the main
argument, the [technical appendices](../output/pdf/focused_appendices.pdf) for
proofs, or the [comprehensive report](../output/pdf/research_report.pdf) for the
full research program. The [teaching guide](../freshman_guide.html) opens locally
in a browser after downloading the package.

This map selects entry points rather than listing every exploratory program.
“Recorded evidence” means a saved research result, not a claim that it was freshly
rerun for publication. [PACKAGE_CHECKS.json](../PACKAGE_CHECKS.json) identifies
the bounded checks actually rerun. An exact computational check can still depend
on an analytic argument; [formal/README.md](../formal/README.md) states the
separate Lean coverage.

| Topic | Question it addresses |
| --- | --- |
| [1. Counterexamples](#1-two-distinct-four-state-counterexamples) | Which witness is being checked or plotted? |
| [2. Controlled onset](#2-controlled-onset-of-cooling-reversal) | How does a reversal appear as a parameter changes? |
| [3. Sharp quarter bound](#3-sharp-four-state-amplitude-bound) | Why is the four-state amplitude supremum 25%? |
| [4. Classification algorithms](#4-asymptotic-classification-and-exact-decisions) | Which asymptotic effects exist for a given generator? |
| [5. Repeated slow modes](#5-degenerate-and-repeated-slow-modes) | What changes when the slow eigenvalue is repeated? |
| [6. Constraints and realizability](#6-rate-constraints-energy-thresholds-and-realizability) | Which modes and effects are compatible with the rate constraints? |
| [7. Higher dimensions](#7-higher-dimensional-constructions-and-the-open-bound) | What survives beyond four states, and what remains open? |
| [8. Finite-time chi-square](#8-finite-time-weighted-euclidean-and-chi-square-ordering) | Can transient reversals occur without eventual reversal? |
| [9. Finite-time entropy](#9-finite-time-relative-entropy-results-and-limits) | Which KL or convex-divergence orderings are proved? |
| [10. Physical-model study](#10-separate-arrhenius-model-and-continuous-robustness-study) | What is certified in the constrained Arrhenius model? |
| [11. Finite strict witness](#11-later-finite-strict-witness-above-8) | Where is the later concrete amplitude witness above 8%? |

## 1. Two distinct four-state counterexamples

Both witnesses invalidate the stated ordinary-cooling triplet implication, but
their generators and finite-time trajectories differ. Do not substitute one
witness's plotted values or certificates for the other.

| Witness | Read and check | Recorded evidence and formal coverage |
| --- | --- | --- |
| Primary witness in the comprehensive research account | [Standalone exact checker](../research/standalone_primary_check.py), [adversarial tests](../research/test_standalone_primary_check.py), [original certificate checker](../research/adversarial_certificate.py) | [Original certificate](../research/adversarial_certificate.json). This is not the Lean companion generator. |
| Rational companion, historically Appendix B and now the focused paper's worked example | [Focused paper source](../focused_paper.tex), [independent arithmetic checker](../research/check_freshman_walkthrough.py), [rational-mode checker](../research/adversarial_rational_mode_certificate.py) | [Arithmetic receipt](../research/detailed_audit/freshman_arithmetic.json), [rational-mode certificate](../research/adversarial_rational_mode_certificate.json), [Lean counterexample](../formal/Mpemba/Counterexample.lean) |

The [focused overview figure](../output/figures/focused_overview.pdf) uses the
rational companion in its relaxation panel. The [Lean relaxation bridge](../formal/Mpemba/Relaxation.lean)
establishes eventual ordering for the specified weighted distance, not a KL
crossing-time certificate.

## 2. Controlled onset of cooling reversal

A continuous family extending the rational companion has an exactly specified
onset threshold, strict triplet failures, and a uniform spectral gap. This is a
classification within that family; it is not a classification of every generator.

- Read: [predictive mechanism](../research/predictive_mechanism.tex).
- Check: [mechanism certificate](../research/mechanism_certificate.py); inspect the
  [recorded certificate](../research/mechanism_certificate.json).
- View: [onset figure](../output/figures/mechanism_onset.pdf).

## 3. Sharp four-state amplitude bound

Under the full triplet-failure premise and a simple slow mode, every finite
cooling pair has amplitude advantage below 25%, and admissible families approach
that supremum. This is an asymptotic amplitude result, not a universal reduction
in relaxation time.

- Read: [full-premise bound](../research/full_triplet_amplitude_bound.tex),
  [quarter-limit family](../research/quarter_limit_family.tex), and
  [single-endpoint criterion](../research/full_triplet_endpoint_criterion.tex).
- Check: [bound identities](../research/full_triplet_amplitude_certificate.py),
  [family certificate](../research/quarter_limit_certificate.py), and
  [thermal witness](../research/quarter_limit_thermal_certificate.py).
- Inspect: [bound receipt](../research/full_triplet_amplitude_certificate.json),
  [family receipt](../research/quarter_limit_certificate.json), and
  [quarter-limit figure](../output/figures/quarter_limit.pdf).
- Formal entry points: [full-generator bridge](../formal/Mpemba/SlowModeBridge.lean)
  and [sharp supremum](../formal/Mpemba/SharpQuarter.lean). The manuscript's stronger
  uniform-gap and all-strict family claims extend beyond that formal coverage.

## 4. Asymptotic classification and exact decisions

Choose an implementation by its input domain. The general simple-mode analysis,
the commensurate-energy exact algorithm, and the full-premise rational algorithm
have different assumptions. They do not classify all finite-time KL or TV crossings.

| Route | Source and canonical implementation | Example evidence |
| --- | --- | --- |
| General four-state simple slow mode | [Analytic classification](../research/universal_classification.tex), [floating classifier](../research/four_state_classifier.py) | [Numerical audit](../research/classifier_audit.json); floating calculations do not certify arbitrary equality cases |
| Rational inputs with commensurate energies, including slow multiplicities | [Commensurate classification](../research/commensurate_classification.tex), [exact classifier](../research/exact_commensurate_classifier.py) | [Example input](../research/commensurate_example.json), [classified output](../research/commensurate_example_classification.json), [regression tests](../research/test_exact_commensurate_classifier.py) |
| Rational rates and bath weights under the full triplet-failure premise, without requiring commensurate energies | [Exact logarithmic endpoint decision](../research/rational_log_endpoint.tex), [exact triplet classifier](../research/exact_rational_triplet_classifier.py) | [Exact-equality fixture](../research/rational_triplet_compact_equality.json), [independent audit](../research/independent_rational_classifier_audit.json) |

The last route has an explicit transcendence-theorem dependency and rejects
inputs outside its premise. Optional resource limits can return “unresolved.”

## 5. Degenerate and repeated slow modes

General repeated-mode asymptotic classification is stated for specified
distances. Under full triplet failure, the separate one-fast-mode theorem gives
a reset-plus-adjacent-pair classification and normal ordering for all convex
divergences at every time. It covers repeated slow multiplicities in four states,
not every repeated slow multiplicity in higher dimensions.

- Read: [degenerate extension](../research/degenerate_extension.tex),
  [four-state full-premise boundary](../research/full_triplet_degenerate.tex), and
  [one-fast-mode classification](../research/one_fast_mode_classification.tex).
- Check: [degenerate identities](../research/degenerate_exact_certificate.py) and
  [one-fast-mode certificate](../research/one_fast_mode_certificate.py).
- Inspect: [degenerate receipt](../research/degenerate_exact_certificate.json),
  [one-fast-mode receipt](../research/one_fast_mode_certificate.json), and the
  [five-state limitation](../research/five_state_degenerate_nonmonotone.py).
  Nonmonotonicity in that limitation is not itself a thermal witness.

## 6. Rate constraints, energy thresholds, and realizability

The weaker G constraints exclude strong effects and asymptotic heating in any
dimension. Energy-span absence theorems and inverse mode-realizability results
have their own premises; G-only sharpness does not establish sharpness under all
F and G constraints.

- Read: [ordered-conductance theorem](../research/ordered_conductance_theorem.tex),
  [energy threshold](../research/energy_threshold.tex), and
  [mode realizability](../research/mode_realizability.tex).
- Check: [conductance identities](../research/ordered_conductance_certificate.py),
  [threshold certificate](../research/energy_threshold_certificate.py), and
  [exact inverse feasibility](../research/exact_mode_realizability.py).
- Inspect: [conductance receipt](../research/ordered_conductance_certificate.json),
  [threshold receipt](../research/energy_threshold_certificate.json), and
  [realizability receipt](../research/exact_mode_realizability.json).

## 7. Higher-dimensional constructions and the open bound

Strict examples extend to every dimension at least four. A separate construction
approaches amplitude advantage `1 − 3/n`; an exact finite five-state witness
exceeds 37.67%. Large energy separations limit practical interpretation. For
`n > 4`, the gap to the upper bound `1 − 2/n` remains unresolved.

- Read: [embedding extension](../research/higher_dimension_extension.tex),
  [amplitude construction](../research/high_dimension_amplitude.tex), and
  [conditional bound and open reduction](../research/high_dimension_bound_reduction.tex).
- Check: [five-state embedding](../research/five_state_extension_certificate.py),
  [clone identities](../research/high_dimension_amplitude_certificate.py), and
  [finite thermal witness](../research/five_state_large_advantage_certificate.py).
- Inspect: [embedding receipt](../research/five_state_extension_certificate.json)
  and [finite-witness receipt](../research/five_state_large_advantage_certificate.json).
  The all-dimensional conclusions depend on the written construction, not finite sampling.
- For the unresolved direction: [graph-coupling LP reduction](../research/general_n_arithmetic_lp.tex)
  and [exact obstruction to one proposed proof route](../research/reset_trace_obstruction.py).
  An obstruction to a proof route does not resolve the open bound.

## 8. Finite-time weighted-Euclidean and chi-square ordering

For a specified thermal pair in a reversible four-state system, the weighted
Euclidean test reduces the finite-time question to endpoints and at most one
interior time. A separate exact example has two crossings while its initial and
eventual ordering remain normal. That example is not under full triplet failure.

- Read: [finite-time classification](../research/finite_chi_classification.tex).
- Check: [transient-only construction](../research/transient_only_chi_certificate.py)
  and the [primary-witness finite-time certificate](../research/finite_time_certificate.py).
- Inspect: [transient receipt](../research/transient_only_chi_certificate.json) and
  [primary finite-time receipt](../research/finite_time_certificate.json).
  These results do not supply a general finite-time KL or TV classification.

## 9. Finite-time relative-entropy results and limits

The near-bath theorem gives a sufficient KL-ordering window under the G
constraints. A different theorem gives all-time heating order for every convex
divergence when the contrast generator is tridiagonal and Metzler. Neither covers
every admissible generator at every temperature.

- Read: [near-bath entropy window](../research/near_bath_entropy_window.tex),
  [tridiagonal contrast theorem](../research/tridiagonal_contrast_entropy.tex), and
  [thermal-ratio obstruction](../research/thermal_ratio_obstruction.tex).
- Check supporting identities: [entropy reductions](../research/check_entropy_reductions.py)
  and [thermal-ratio certificate](../research/thermal_ratio_obstruction_certificate.py).
- Inspect: [reduction receipt](../research/detailed_audit/entropy_reduction_identities.json)
  and [thermal-ratio receipt](../research/thermal_ratio_obstruction_certificate.json).
  The identity checks support the analytic arguments; they do not independently
  prove the universal theorems. The obstruction is not a demonstrated KL reversal.

## 10. Separate Arrhenius-model and continuous robustness study

This model has continuous certificates for the specified KL-distance advantage
and threshold-time savings throughout its simultaneous ±5% parameter box. Two
triplets pass, so this is separate from the all-triplet-failure counterexample.
The thermal-collision realization is an idealized stochastic model, not material
or experimental validation.

- Read: [model and finite-time results](../research/physical_finite_time_results.tex),
  [robustness proof and reductions](../research/physical_robustness_certificate.tex),
  and [thermal-collision realization](../research/thermal_collision_realization.tex).
- Reproduce: [continuous certificate runner](../research/run_robust_finite_time_certificate.py)
  and [thermal-collision model](../research/thermal_collision_model.py).
  The continuous runner requires C++17 and MPFR/GMP in addition to Python.
- Inspect: [eight-obligation certificate](../research/robust_finite_time_certificate/certificate.json),
  [collision audit](../research/thermal_collision_audit.json), and
  [finite-time figure](../output/figures/physical_finite_time.pdf).
- Distinguish the [activated-intermediate study](../research/activated_intermediate_study.tex):
  it is a numerical coarse-graining investigation, not a continuous robustness
  certificate. Earlier incomplete enclosure attempts remain historical evidence.

## 11. Later finite strict witness above 8%

The follow-on study uses a finite member of the quarter-limit family with a
strict perturbation. Its recorded amplitude enclosure is approximately
8.1718–8.1823%, using exact arithmetic and a stated spectral perturbation bound.
The associated thermal screening is high-precision numerical evidence. This is
not a hardware implementation or a physical tolerance-box certificate.

- Check: [finite strict certificate](../research/home_counterexample/certify_intermediate.py).
- Inspect: [certificate receipt](../research/home_counterexample/intermediate_certificate.json)
  and the [fresh bounded replay result](../PACKAGE_CHECKS.json), which records an
  identical regenerated receipt.
- Explore separately: [numerical screening script](../research/home_counterexample/check_intermediate.py)
  and [screening results](../research/home_counterexample/intermediate_results.json).

[Back to project home](../README.md) · [Documentation hub](README.md) · [Detailed command reference](../RUNME.txt)
