# Independent skeptical audit of the certified counterexample

Date: 2026-10-11. Verdict: the final rational-mode example disproves the stated four-state inference for asymptotic cooling at a finite positive bath temperature. This verdict follows from exact algebra and certified inequalities, not agreement between numerical searches.

## Independently reconstructed data

Set beta_b=1, q=(1,199/200,2/5,1/400), and E_i=-log(q_i). These are four distinct, strictly increasing energies. Let the six downward rates, in order (12,13,14,23,24,34), be

\[
\left(1,28,\frac{249}{10},\frac{330154675693}{1470375000},
\frac{4983}{200},\frac{1527327066}{151990625}\right).
\]

For i<j set M_ij to this downward rate and M_ji=M_ij q_j/q_i; set each diagonal to minus its off-diagonal column sum. Independent reconstruction used no discovery-code imports.

The exact slow eigenfunction is

\[
w=(1,-759/1000,-48957/80000,-1/250)^T,
\qquad\lambda_2=-4966671/250000.
\]

Verified exactly: every off-diagonal rate is positive, columns sum to zero, detailed balance holds, M^T w=lambda_2 w, and q^T w=0. With

\[
N=\sum_iq_iw_i^2=27567998009/16000000000,
\qquad v=\operatorname{diag}(q)w/N,
\]

Mv=lambda_2 v and w^T v=1 hold exactly. Thus w^T pi(beta) is precisely the coefficient in this chosen normalization.

## Spectral and predicate checks

The independently computed characteristic polynomial is

\[
\frac{x(250000x+4966671)}{632507162050000000000000}
\left(2530028648200000000x^2+
1002356769715761942000x+50981800562345565497767\right).
\]

Exact Sturm counts locate one root of the quadratic in (-337,-336) and one in (-60,-59). Since the rational root is -19.866684, it is the unique slowest nonzero root, separated by a substantial gap. The independent script evaluates the original branch-dependent predicate, not merely its simplified equivalent. All four predicates fail, with F<0 and G<0 for each triplet; all inequalities are strict.

Reproduction: `python3 research/independent_mode_audit.py`; executed output is saved in `research/independent_mode_audit.log`.

## Thermal witness and rigorous interval audit

The separately authored `research/adversarial_rational_mode_certificate.py` was read line by line. Its interval arithmetic uses exact rationals. The floating-point approximation only proposes rational endpoints for q_i^(1/n); acceptance verifies lower^n <= q_i <= upper^n exactly, with nonnegative endpoints. This correctly encloses the unique positive roots. Interval addition, multiplication, and reciprocal operations are valid, and the denominator interval is positive.

The two inverse temperatures are beta_h=1/1000 and beta_c=1/100, strictly between zero and beta_b. The saved exact rational bounds imply the following deliberately rounded-outward, terminating-decimal rational enclosures:

\[
-0.093755177882 < b(1/1000) < -0.093755177881,
\]

\[
-0.093857287853 < b(1/100) < -0.093857287852.
\]

Both overlaps are negative. Hence |b(1/1000)|<|b(1/100)|; the exact certificate's margin lower bound exceeds 0.0001021099703224. This conversion from a signed interval to the absolute-overlap inequality was checked. The temperatures are finite (T_h=1000 and T_c=100 in units k_B=1), and T_b=1.

With the simple slow mode already verified, spectral expansion proves eventual reversal for every fixed norm and for relative entropy. This is an ordinary asymptotic Mpemba effect: neither displayed overlap vanishes. No unproved claim about a turning point or local maximum is needed. The source's three-state tangent construction is only the derivation of its necessary predicate; the four-state definition is directly met by this finite thermal pair.

## Robustness and scope

The positive-rate, strict-predicate, simple-eigenvalue and strict-overlap inequalities are all open conditions. Continuity of the reversible simple eigenprojection therefore supplies an open neighborhood of counterexamples in the unrestricted four-state parameter space. This is not an isolated equality case, a zero-rate limit, an infinite bath-temperature limit, or a spectral-degeneracy artifact.

The finding is a counterexample to the exact full-four-state implication stated in the assignment and in the source's conjecture. It does not undermine the source's explicit zero-temperature strong-violation theorem. Nor does it imply that any three-state necessary condition is false: combining triplets into a four-state relaxation changes the global slow eigenfunction and thermal overlap.
