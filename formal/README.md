# Lean companion to Appendix B

The exact rational-mode counterexample and a connection to actual relaxation
are formalized in Lean 4.19.0 with pinned mathlib commit
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. The sharp quarter-amplitude
supremum is also formalized, including a valid limiting family. No unfinished
proof is admitted in the built library.

Run from this directory:

```sh
./check.sh
```

The script uses the workspace-local toolchain when available, otherwise an
installed `lake`. It builds every module, prints the axiom dependencies of
43 named theorems, and rejects any axiom beyond `propext`, `Classical.choice`, and
`Quot.sound`. The mathlib dependency and its dependencies are fixed by
`lake-manifest.json`; local toolchains and downloaded build files are ignored.

## Verified counterexample

`Mpemba.appendixB_counterexample` in `Mpemba/Counterexample.lean` collects:

* Positive ordered bath weights and distinct ordered energies `energy i = -log(q i)`.
* All twelve positive transition rates, column conservation, and detailed balance.
* Failure of all four original triplet predicates, with the original inclusive
  kappa branches and strict sign tests.
* The exact left and right slow eigenvectors and their normalization.
* The characteristic polynomial, algebraic multiplicity one of the rational
  slow eigenvalue, and dominance over every other nonstationary real root.
* Finite cooling preparations `0 < 1/1000 < 1/100 < 1` and strict Gibbs-overlap reversal.
* An initial-to-eventual reversal under the printed generator's master equation.

The real Gibbs weights are proved to equal `exp(-β * energy i)`. Thermal
enclosures are established by rational integer-power inequalities in Lean;
no decimal approximation or external certificate is accepted as a premise.
The formally checked overlap intervals are deliberately wider than the
manuscript's printed numerical intervals:

```text
-0.09376 < b(1/1000) < -0.09375
-0.09386 < b(1/100)  < -0.09385
```

`Data.lean` contains the literal Appendix B data; states 1–4 use indices 0–3.
`M i j` is a transition from j to i. `Multiplicity.lean` establishes algebraic
simplicity, in addition to the eigenspace-line statement in `Spectral.lean`.

## Physical bridge and its scope

`Dynamics.lean` constructs the evolution from the stationary, slow, and two
fast components. Lean verifies the initial condition and the master equation
coordinate by coordinate. `Uniqueness.lean` proves that this is the unique
global differentiable solution with that initial condition.

The distance used in the current formal bridge is

```text
D(p) = Σ_i (p_i - π_bi)^2 / q_i,   π_bi = q_i / Q,   Q = Σ_i q_i.
```

Multiplication by Q gives Pearson chi-square divergence. This also orders
the associated weighted Euclidean norm. The bridge is proved for this
specified distance; a relative-entropy/KL bridge is not formalized here.

`Relaxation.lean` proves for arbitrary normalized initial vectors that strict
slow-overlap ordering gives strict distance ordering at all sufficiently
late times. The spectral expansion and its limiting estimates are proved,
not supplied as assumptions. `InitialDistance.lean` separately proves

```text
D(π(1/100)) < 24 < D(π(1/1000)).
```

Thus the hotter state starts farther away and eventually becomes closer.
`master_equation_distance_reversal` transfers this conclusion to any global
solutions with the two specified Gibbs initial conditions.
`reversal_at_a_finite_positive_time` establishes existence at a finite
positive real time; no numerical crossing time is certified.

## Sharp quarter-amplitude supremum

`QuarterGibbs.lean` proves the analytic bound from explicit mode-shape and
overlap-sign hypotheses. It includes monotonicity of the conditional
three-state Gibbs mean and bounds on the first-three probability mass.
`QuarterBound.lean` proves the resulting ratio below 4/3 and amplitude
advantage below 1/4.

`TripletLogic.lean` proves the full original-branch equivalence to the F,G
cone. `QuarterBoundary.lean` checks the three reset-boundary polynomial sign
cases, including zero-rate faces. `QuarterSchur.lean` checks the Schur
numerator identity, its sign reduction, and the eigenvector ratio comparison.

`FourState.lean` defines the general four-state reversible generator at
ordered bath weights `(1,t,u,v)`, its weighted Helmert basis, and its contrast
matrix. `TopMode.lean` proves the nonnegative orientation of a simple top
mode from a positive-semidefinite variational certificate. `Reset.lean`
proves reset subtraction preserves the F/G cone and eigenvectors, leaves
nonnegative rates with a zero first-row rate, and gives a nonpositive shifted
eigenvalue by a Dirichlet-form identity. `SchurSpectral.lean` proves the strict
diagonal gaps needed for elimination. `ModeShape.lean` combines these facts
with the boundary polynomial to derive all three required shape inequalities.

`OriginalTriplets.lean` proves that the four original predicates, including
the inclusive kappa branches, are equivalent to the eight denominator-cleared
F/G inequalities. `CoolingOverlap.lean` derives strict negative cooling
overlaps from nonnegative contrast coordinates and removes the orientation
and normalization choices using simplicity.

The full-generator entry points in `SlowModeBridge.lean` are:

* `full_generator_quarter_bound`: the cold/hot absolute overlap ratio is
  strictly less than `4/3`.
* `full_generator_amplitude_advantage_below_quarter`: `1 - hot/cold`
  absolute amplitude advantage is strictly less than `1/4`.

These theorems assume positive downward rates, `0<v<u<t<1`, failure of all
four original predicates, and `SimpleSlowMode` for the full generator. That
structure explicitly states a nonzero bath-orthogonal left eigenfunction,
its eigen-equation, one-dimensional eigenspace, and the usual variational
largest-eigenvalue condition on the bath-orthogonal subspace. The bridge
proves the contrast-coordinate representation and transfers those conditions;
mode shape and overlap signs are conclusions, not assumptions. The bath
inverse temperature is normalized to one and preparations satisfy
`0 < hot ≤ cold < 1`.

`SharpQuarter.lean` defines `attainableAdvantages` from positive rates,
ordered bath weights, all four original triplet failures, a simple slow mode,
finite strictly ordered cooling preparations, and strict ordinary cooling
amplitude reversal. Its main results are:

* `sharp_quarter_isLUB`: `1/4` is the least upper bound of this set.
* `sharp_quarter_not_attained`: no member equals `1/4`.
* `exists_advantage_arbitrarily_close`: for every positive tolerance there is
  an admissible advantage strictly within that tolerance below `1/4`.

The construction in `SharpFamily.lean` uses a simpler rational family than
the manuscript's uniformly gapped family. It retains

```text
q = (1, 1-epsilon/3, 3/5, epsilon^2)
w = (1, -4/5, -1/3-4*epsilon/9+5*epsilon^3/3, -epsilon)
lambda = -3
```

It fixes the downward rate `a12=1` and gives explicit rational functions for
the remaining five rates. Lean verifies bath orthogonality and the exact
eigen-equation. Thirteen strict checks at epsilon zero, together with proved
continuity, give a nonempty interval of positive epsilon on which all rates
are positive and every original triplet fails. Two cleared inequalities,
F124 and G234, are exact equalities. This is allowed by the original strict
predicate tests. `PositiveMode.lean` proves that the positive contrast mode
and connecting positive edges certify the simple largest contrast eigenvalue;
`SlowModeBridge.lean` transfers it to the full generator.

`SharpThermal.lean` sets `epsilon=exp(-L^2)`, `hot=L^(-3)`, `cold=L^(-1)`
for sufficiently large finite L. It checks `0<hot<cold<1` and proves

```text
b_hot  -> -1/30
b_cold -> -2/45
1 - |b_hot|/|b_cold| -> 1/4
```

`eventually_family_attainable` connects these limits to the full admissible
generator conditions. The sharp supremum therefore concerns actual valid
members of the family, not only a conditional overlap formula.
`FourStatePhysics.lean` additionally verifies conservation, detailed balance,
all twelve positive off-diagonal rates, and the real-exponential Gibbs
interpretation for the general generator parametrization.

This verifies the sharp amplitude supremum. The manuscript's stronger uniform
gap estimate, all-eight-strict perturbed sharpness family, general KL bridge,
and time-advance corollaries are outside the formal coverage. The 25% quantity
is an asymptotic amplitude advantage, not a claimed percentage of cooling time.

## First setup on another machine

Install elan, Git, Python 3 and a POSIX shell. Keep `lean-toolchain`,
`lakefile.toml` and `lake-manifest.json` with the sources. From the package root:

```sh
cd formal
lake exe cache get
./check.sh
```

Elan selects Lean 4.19.0; Lake fetches the dependency revisions in the committed
manifest, and the cache command downloads precompiled mathlib dependencies.
Initial setup requires network access. Do not run `lake update` when reproducing
the committed dependency set. Once dependencies are available, `./check.sh`
is the verification command.

On this Mac, the compiled mathlib cache downloader encountered a dyld linker
error. Running its Lean interpreter entry point worked:

```sh
lake env lean --run .lake/packages/mathlib/Cache/Main.lean get \
  Mathlib.LinearAlgebra.Matrix.Charpoly.Coeff \
  Mathlib.LinearAlgebra.Matrix.PosDef \
  Mathlib.Algebra.Polynomial.FieldDivision \
  Mathlib.Analysis.Normed.Module.FiniteDimension \
  Mathlib.Analysis.ODE.Gronwall \
  Mathlib.Analysis.Calculus.Deriv.Prod \
  Mathlib.Analysis.SpecialFunctions.Pow.Real \
  Mathlib.Analysis.SpecialFunctions.Pow.Continuity \
  Mathlib.Analysis.SpecialFunctions.ExpDeriv \
  Mathlib.Data.Matrix.Notation Mathlib.Tactic.FinCases \
  Mathlib.Tactic.NormNum Mathlib.Tactic.Linarith \
  Mathlib.Tactic.FieldSimp Mathlib.Tactic.Ring Mathlib.Tactic.Positivity \
  Mathlib.Tactic.FunProp
```

The proofs use standard kernel-checked tactics. They contain no `sorry`,
custom axioms, or `native_decide` dependencies. Formal correctness does not
settle literature priority or the correspondence of the published conjecture
to its transcription; those still require scholarly review.
