# Source and definition audit

Checked 2026-10-11. This is a primary-source audit and independent algebraic reconstruction, not a claim of an exhaustive literature review.

## Sources and status

1. Avitan, Factor and Gelbwaser-Klimovsky, [arXiv:2603.04567](https://arxiv.org/abs/2603.04567), v1 submitted 4 March 2026, 19:54:17 UTC. The abstract page listed only v1. [Full text and supporting information](https://arxiv.org/html/2603.04567v1) were read, especially Eq. (3), the N-level protocol, and SI Eqs. S2–S14. The finite-temperature N-level inference is explicitly a conjecture supported by simulations. The analytic theorem concerns zero-temperature baths, nondegenerate generator eigenvalues, and the stronger hypothesis C>A+B on every triplet. That theorem does not settle the requested finite-temperature premise. The numerical discussion includes equal and zero rates; it supplies no general boundary proof.
2. Gelbwaser-Klimovsky, [arXiv:2609.09988](https://arxiv.org/abs/2609.09988), v1 submitted 9 September 2026, 10:15:30 UTC; only v1 was listed. [Full text and SI](https://arxiv.org/html/2609.09988v1) derive exact three-state conditions and propose a search near a fully degenerate relaxation generator. Section II repeats Eq. (3) and attributes the N-level inference to the earlier work; it adds no general proof. Its approximately 2% four-state effects with no exact three-state mechanism do **not** disprove the earlier predicate, which is weaker than the exact three-state condition. Its sampling used 500,000 realizations, rates from 0.001 to 1, and dimensionless energies from 0 to 5.
3. Lu and Raz, [PNAS 114, 5083–5088 (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5441807/), full text read. Their definition requires eventual distance reversal persisting for all sufficiently large times, and their simple-slow-mode overlap criterion is sufficient. They treat direct cooling and inverse heating. They require distances to decrease during relaxation, increase along the initial Gibbs family away from the bath, and be continuous and convex. They name relative entropy and L1 distance as examples. The scalar criterion explicitly assumes a spectral gap between the first two nonzero modes.
4. Klich, Raz, Hirschberg and Vucelja, [Phys. Rev. X 9, 021060 (2019)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.021060), abstract and publication record checked; [arXiv:1711.05829](https://arxiv.org/abs/1711.05829). This introduces the strong effect: specially chosen initial temperatures eliminate the slow relaxation contribution and produce an exponentially faster approach.

Searches for the original identifier, title, conjecture, counterexample and corrections found no additional primary-source resolution. The September paper was inspected rather than treated as a resolution from its abstract. This bounds the checked evidence; it cannot exclude an unindexed or unpublished result. The mathematical target remains unchanged.

## Exact triplet predicate and independent algebra

With the ordered triplet i<j<k, the source's symbols (-,0,+) correspond to (i,j,k). Its rates are destination-first, consistent with column probability vectors. The goal's A=a_ik, B=a_jk, C=a_ij, t=exp[-beta_b(E_j-E_i)], u=exp[-beta_b(E_k-E_j)] give exactly the two branches printed in Eq. (3) and S9:

\[
\kappa=(A+B)/(C(1+t)),\quad
F=A(1+tu/2)-B(1+u/2)-C(1-t),\quad G=At-B.
\]

The proposed predicate is the faithful transcription. The strict inequalities F>0 and G>0 matter. SI S11–S14 explicitly use the original off-diagonal rates of the full four-state matrix. A triplet's diagonal cannot simply be copied from a four-state principal submatrix if one means an isolated conservative generator.

There is a useful exact simplification of the branch logic. Define

\[
H=A+B-C(1+t),\quad d=(1-t)/(1+t)>0,\quad h=2/(1+t)+u/2>0.
\]

Direct expansion gives

\[
F=hG+dH.
\]

For H>0, failure means F<=0, which forces G<0. For H<0, failure means G<=0, which forces F<0. At H=0, F=hG, so the two branches agree. Consequently, throughout the stated domain,

\[
\boxed{\neg P_{ijk}\iff F\le0\ \text{and}\ G\le0.}
\]

This was independently checked by expanding the A, B and C coefficients. Thus, at fixed energies, all four failures form the intersection of eight closed linear half-spaces with the positive rate orthant: one convex polyhedral cone, rather than a genuine union of branch cones.

The SI's fast-eigenvector formulas S4/S8 divide by a_+0-a_+-=-uG and do not explicitly cover its zero. The predicate itself is well-defined there when C>0. If G=0, F=dH. If additionally H=0, A=C and B=At; all nonzero eigenvalues of the isolated triplet coincide. If G=0 but H differs from zero, directly diagonalizing the conservative triplet gives its two nonzero eigenvalues

\[
-A(1+t+tu),\qquad -[C(1+t)+Atu].
\]

Thus that denominator-zero set is not automatically a degenerate spectrum, and G=0 alone does not imply failure of P. No division by G is needed for the boxed equivalence.

At kappa=1, equality F=G=0 is precisely the degenerate triplet just described. At kappa=1 with G<0 both tests fail; with G>0 both pass. There is no contradictory branch assignment at this boundary.

## Scope of the target and definition repairs

The original article's discussion of a maximum followed by “i.e. derivative equals zero” should not be used as a general equivalence. A stationary point can fail to be a local maximum, and a stationary point at the bath is excluded. In four states the nature of a turning point requires proof. The September SI explicitly addresses the stationary-point issue in its three-state setting; that argument must not silently be generalized.

Use beta>0 for finite positive initial temperatures. Cooling means 0<beta_h<beta_c<beta_b; heating means beta_b<beta_h<beta_c<infinity. Cooling reversal is |b(beta_h)|<|b(beta_c)|. Heating reversal is |b(beta_c)|<|b(beta_h)|. Test the two bath sides separately. Infinite-temperature and zero-temperature preparations can be useful limits but are not necessary witnesses for this task.

For a simple nonzero slow eigenvalue, let Mv=lambda v and w^T M=lambda w^T, normalized by w^T v=1. Since M pi_b=0 and lambda differs from zero,

\[
0=w^TM\pi_b=\lambda w^T\pi_b,
\]

so b(beta)=w^T[pi(beta)-pi_b]=w^T pi(beta) exactly. Any constant rescaling of w multiplies all overlaps equally and does not change their absolute ordering.

Independently of distance conventions, a strict absolute-overlap inequality certifies eventual reversal for every fixed norm:

\[
p(t;\beta)-\pi_b=e^{\lambda_2t}b(\beta)v_2+O(e^{\lambda_3t}).
\]

It also certifies eventual relative-entropy reversal, because its leading nonzero term is

\[
D_{KL}(p(t;\beta)\Vert\pi_b)
=\tfrac12 e^{2\lambda_2t}b(\beta)^2\sum_i v_{2i}^2/\pi_{bi}
+o(e^{2\lambda_2t}).
\]

If the farther preparation has b=0 while the nearer preparation does not, the faster asymptotic exponential decay is the strong effect. These statements do not assert a prescribed finite crossing time, nor exclude a transient crossing when asymptotic ordering is ordinary.

For an arbitrary nonsymmetric distance, absolute coefficient ordering alone needs additional hypotheses; the explicit norm and relative-entropy arguments avoid that subtlety. Relative entropy of initial Gibbs distributions is strictly increasing away from the bath on each side, since

\[
\frac{d}{d\beta}D_{KL}(\pi(\beta)\Vert\pi_b)
=(\beta-\beta_b)\operatorname{Var}_{\pi(\beta)}(E).
\]

With a degenerate slow eigenvalue, replace b v by the full spectral projection onto the slow eigenspace; a chosen basis vector's coefficient is not intrinsic. Norm/entropy comparisons can then depend on the entire projected vector. Detailed balance still guarantees diagonalizability, but not scalar simplicity. Zero rates require connectedness/irreducibility for a unique bath equilibrium, and C=0 makes the printed kappa undefined. One may study a carefully justified closure via F,G, but that is an extension of the printed quotient test, not a literal evaluation at C=0.

## Additional source caveat

The original main text calls (i) C>A+B and (ii) A+B>C>A-B two ways of failing the zero-temperature test, then says (i) implies (ii). As printed, this implication is false because their comparisons with A+B are opposite. The SI's weaker condition is just C>A-B, which *is* implied by C>A+B. This textual inconsistency does not affect the finite-temperature Eq. (3) transcription or the strong-violation theorem's explicit hypothesis.
