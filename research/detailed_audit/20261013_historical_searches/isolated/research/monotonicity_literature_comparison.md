# Monotonicity literature comparison

Checked 11 October 2026: [Ohga, Hayakawa and Ito, arXiv:2410.06623v2](https://arxiv.org/abs/2410.06623), submitted 19 October 2024. The abstract history lists v2 as latest. Read the [full text](https://arxiv.org/html/2410.06623v2), including Result 2 and SI B, C and D—not only the abstract. Local source copies are `sources/2410.06623v2.html` and `.txt`. The HTML's body displays a different rendered date; cite the version-history date above.

## Exact comparison

| Issue | OHI Result 2 and SI D | Present report's original order theorem |
|---|---|---|
| State space | Finite, partially ordered | Finite, totally ordered by distinct energies |
| Dynamics | Irreducible, Gibbs-stationary; detailed balance unnecessary | Positive reversible rates |
| Order assumptions | Energy monotone; jumps connect comparable states; monotone coupling | Adjacent-difference matrix K is Metzler |
| Bath hypothesis | Conditions throughout a temperature interval | One fixed positive finite bath temperature |
| Spectral hypothesis | None | Simple slowest nonzero eigenvalue |
| Conclusion | No cooling energy-expectation crossing at any time | No asymptotic cooling/heating reversal in norms or KL |

SI D2 proves monotonicity from coupling (S4), and explicitly treats total orders directly. D3–D4 use Harris association preservation (S5), equilibrium association (S6), and covariance (S37) for partial orders. SI B extends **Result 1**, not Result 2, to KL. Accordingly, do not attribute an all-time KL no-Mpemba theorem to their printed Result 2. [Source](https://arxiv.org/html/2410.06623v2#A4).

## Mathematical relationship: independently checked deductions

These deductions identify how the report fits the established monotonicity framework; they are not a claim of priority.

Let L=M^T, Bf=(f_{i+1}-f_i), and K=BLH_0 as in the report. The identity BL=KB implies

\[
B e^{tL}=e^{tK}B.
\]

Therefore K Metzler implies preservation of increasing functions for every t>=0. Conversely, preservation applied to a function whose adjacent differences are a coordinate vector, and differentiation at zero, implies every off-diagonal K entry is nonnegative. Thus the matrix criterion is **exactly stochastic monotonicity on this finite total order**, not a different mechanism.

A monotone coupling can also be constructed without citing an unproved equivalence. Choose r sufficiently large that P=I+L/r is stochastic and I+K/r is entrywise nonnegative. Then P preserves increasing functions. Its transition rows are stochastically ordered, so common-uniform inverse-CDF updates preserve the order of two starting states. Apply those updates at a common Poisson clock of rate r. Each marginal has generator L, and the coupled process never crosses. This identifies the report's hypotheses with the total-order coupling framework.

On a total order all jumps connect comparable states automatically, and every probability law has nonnegative covariance for two increasing functions. Thus the association machinery needed for a general partial order is unnecessary here. In particular, no bath-dependent family of generators is required: at the fixed bath, e^{tL}E is increasing and

\[
\partial_\beta\bigl[\pi(\beta)^T e^{tL}E\bigr]
=-\operatorname{Cov}_{\pi(\beta)}(e^{tL}E,E)\le0.
\]

This proves ordinary energy ordering at **all times**, directly from K Metzler. It works for heating as well as cooling, and needs neither reversibility nor a spectral simplicity assumption (although Gibbs stationarity identifies the intended equilibrium).

## Stronger distributional consequence, independently checked

The theorem agent's proposed all-time convex-divergence extension is correct under the report's reversible total-order assumptions. It is an elementary consequence of stochastic monotonicity plus reversibility, and should be presented with a proof and without a general novelty claim.

Write p_beta(t)=e^{tM}pi(beta), r_beta(t)=p_beta(t)/pi_b componentwise. Reversibility gives

\[
r_\beta(t)=e^{tL}r_\beta(0),\qquad
r_{\beta,i}(0)\propto e^{(\beta_b-\beta)E_i}.
\]

For cooling the likelihood ratio is increasing; for heating it is decreasing. These properties survive e^{tL}. Gibbs preparations are stochastically ordered by temperature, and that order also survives the semigroup.

For cooling, p_h(t)>=st p_c(t)>=st pi_b. For a convex function Phi finite on positive arguments, select a nondecreasing subgradient g. The supporting-line inequality gives

\[
D_\Phi(p_h\Vert\pi_b)-D_\Phi(p_c\Vert\pi_b)
\ge\sum_i[p_{h,i}-p_{c,i}]g(r_{c,i})\ge0,
\]

where D_Phi(p||pi_b)=sum_i pi_bi Phi(p_i/pi_bi). The second inequality follows from stochastic order and the increasing function g(r_c). For heating, p_c<=st p_h<=st pi_b and r_h is decreasing; applying the supporting line at r_h gives D_Phi(p_c||pi_b)>=D_Phi(p_h||pi_b). Positive initial Gibbs weights and positive pi_b keep ratios inside the domain. No slow-eigenvalue simplicity is needed.

As a second check, for differentiable Phi,

\[
\partial_\beta D_\Phi(p_\beta(t)\Vert\pi_b)
=-\operatorname{Cov}_{\pi(\beta)}
\bigl(e^{tL}\Phi'(r_\beta(t)),E\bigr).
\]

It is nonpositive for cooling and nonnegative for heating. Taking Phi(z)=z log z proves all-time KL ordering. Taking Phi(z)=|z-1| proves all-time L1/total-variation ordering via the subgradient argument. **This does not prove all-time ordering for every arbitrary norm**; the report's spectral argument supplies the asymptotic every-norm statement under its additional simple-mode assumption.

## Broader novelty assessment (11 October 2026)

The report now includes an integrated broader comparison in section 4.5, generated
from `research/novelty_assessment.tex`. This supersedes the earlier single-paper
scope and proposed wording. The report retains the distinctions and proofs above.

- **Fill and Kahn**, arXiv:1109.6075v3 (6 September 2013), sections 3.1–3.3:
  Proposition 3.2, Corollary 3.3 and equation (3.2) provide comparison and
  density-ratio precedents. Their different-kernel/common-initial-law setup is
  distinguished from our one-kernel/different-Gibbs-preparation setup. Their
  printed corollary is not attributed the exact all-convex result used here.
- **Baake and Sumner**, arXiv:2007.11433v3 (25 September 2022), section 2.1,
  Fact 3.1 and Proposition 3.11: direct precedents for equal-input reset generators
  and the transformed-generator characterization. The report shows explicitly
  how their cumulative-coordinate matrix has our K as its lower block.
- **Vu and Hayakawa**, arXiv:2410.06686v2 (6 January 2025), main-text Lorenz
  construction, conditions (a)–(c), Results 2–4: existing thermomajorization gives
  the all-convex-divergence characterization. The report derives the common-knot
  Lorenz comparison from its stochastic and likelihood-ratio ordering hypotheses.
  This is a standard-method application, not an independently novel principle.

The counterexample is the central contribution. The all-lower-kappa and
nonincreasing gap-law statements remain proved applications; their individual
priority is not established. The cone identity and reset compatibility remain
useful problem-specific deductions without separate priority claims. No claim
of a universal-monotone reversal is added to the counterexample.

Versioned local source copies and exact inspection scopes are recorded in
`research/novelty_search_audit.json`. The review is broader but not exhaustive;
not locating an identical theorem is not proof of originality.
