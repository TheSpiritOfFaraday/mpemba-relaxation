# Historical proof-oriented research notes

These notes retain an early stage of the investigation. Their continuation questions describe that stage; consult the comprehensive report and later certificates for the current results and remaining open problems.

Date: 2026-10-11. Scope: finite positive bath inverse temperature, strictly ordered energies, positive rates, reversible generator, asymptotic simple slow-mode Mpemba effect. These are derivations performed during this project, not claims of literature novelty.

## 1. The triplet cone and its boundary

Write H=A+B-C(1+t), G=At-B, and F=A(1+tu/2)-B(1+u/2)-C(1-t). Direct algebra gives

F = h G + d H,  h=2/(1+t)+u/2>0, d=(1-t)/(1+t)>0.

For H<0 the paper's failure condition is G<=0, which implies F<0. For H>0 its failure condition is F<=0, which implies G<0. For H=0, F=hG and either branch gives failure exactly when G<=0. Consequently, for distinct energies and finite positive bath temperature,

**NOT P if and only if F<=0 and G<=0.**

Thus the four-triplet feasible region is one convex polyhedral cone in the six downward rates, not a genuinely disjoint union of sixteen cones. The branch language remains useful for identifying structured subfamilies. Eq. (3) and Supporting Information Eq. (S9) in the downloaded v1 text have exactly the proposed transcription. Its >= and <= branches overlap at k=1, where the two strict tests are equivalent.

## 2. Adjacent-difference representation

Let L=M^T be the backward generator acting on functions f. Thus

(L f)_i = sum_{j != i} a_{j i}(f_j-f_i).

Set y_i=f_{i+1}-f_i, for i=1,...,n-1. Define H_0 by (H_0 y)_j=sum_{l<j}y_l, so that f=f_1 1+H_0 y. Let B f=y. Since L1=0,

B L f=K y, where K=B L H_0.

Explicitly,

K_{i,l}=sum_{j>l}(L_{i+1,j}-L_{i,j}).

For off-diagonal entries this is

- l<i: K_{i,l}=sum_{j<=l}(a_{j,i}-a_{j,i+1});
- l>i: K_{i,l}=sum_{j>l}(a_{j,i+1}-a_{j,i}).

The diagonal is

K_{i,i}=-sum_{j<=i}a_{j,i+1}-sum_{j>i}a_{j,i}.

In the basis consisting of the constant function and the columns of H_0, L is block triangular with diagonal blocks 0 and K. Therefore K has exactly the nonzero generator eigenvalues. No equilibrium weighting or numerical diagonalization is needed for this assertion.

For four states its six off-diagonal entries are

K12=a32+a42-a31-a41; K13=a42-a41;
K21=a12-a13; K23=a43-a42;
K31=a13-a14; K32=a13+a23-a14-a24.

Here the a indices retain the original destination,source convention.

## 3. Theorem: a directly checkable finite-temperature sufficient condition

**Theorem.** Consider an n-state reversible generator with strictly positive rates, distinct ordered energies, a finite positive bath inverse temperature, and a simple slowest nonzero eigenvalue. If all off-diagonal entries of K above are nonnegative, then neither an asymptotic cooling nor an asymptotic heating Mpemba effect occurs for positive initial temperatures.

**Proof.** K is Metzler (its off-diagonal entries are nonnegative). The Perron theorem applied after a sufficiently large scalar diagonal shift gives a nonnegative, nonzero right eigenvector y for its spectral abscissa. Because the spectrum of K is precisely the nonzero spectrum of the reversible generator, this spectral abscissa is lambda2. Solve B w=y with the additive constant chosen so pi_b^T w=0. The relation BLw=K Bw=lambda2 Bw implies Lw-lambda2 w is constant. Its bath expectation is zero because pi_b^T L=0 and pi_b^T w=0, so that constant is zero. Thus w is a left slow eigenvector of M, and w is nondecreasing in energy and nonconstant.

For b(beta)=sum_i w_i pi_i(beta), fixed w and fixed bath,

b'(beta)=-Cov_{pi(beta)}(w,E)
=-sum_{i<j} pi_i(beta) pi_j(beta) (w_j-w_i)(E_j-E_i)<0.

Strictness follows because all Gibbs probabilities are positive, energies are strictly ordered, and at least one w difference is positive. Since b(beta_b)=0, b(beta)>0 for beta<beta_b and b(beta)<0 for beta>beta_b. Thus |b| is strictly decreasing on the cooling interval 0<beta<beta_b and strictly increasing on the heating interval beta_b<beta<infinity. These are precisely the normal asymptotic orderings. Multiplying w by a nonzero normalization scalar does not change this conclusion. QED.

Irreducibility of K is not required here: its Perron vector may have zero components, but is nonzero. Simplicity of lambda2 makes this the unique slow mode up to scale. This theorem does not independently assert absence of every possible finite-time crossing for every possible distance functional.

## 4. Corollary: all lower k branches, any number of states

**Corollary.** If every energy-ordered triplet satisfies k<=1 and fails the predicate, the preceding theorem applies.

Indeed H<=0 and G<=0 give

A+B<=C(1+t), B>=At, hence A<=C.

Thus for j<i, triplet (j,i,i+1) gives a_{j,i+1}<=a_{j,i}. For j>i+1, triplet (i,i+1,j) gives a_{j,i}<=a_{j,i+1} by detailed balance and G<=0. Every summand in every off-diagonal K entry is therefore nonnegative.

This establishes the intended absence implication on the entire all-k<=1 chamber at positive finite bath temperature, including its equality faces. It does not assume equally spaced energies or a fixed energy span. It generalizes directly to n states. What prevents extending this proof to all admissible rates is that failure with k>1 does not force A<=C.

## 5. Corollary: decreasing downward gap laws, any number of states

Suppose a_{ij}=f(E_j-E_i) for i<j, with f strictly positive and nonincreasing on positive gaps. Detailed balance fixes upward rates as a_{ji}=exp[-beta_b(E_j-E_i)]f(E_j-E_i). Both the downward function f and upward function exp[-beta_b x] f(x) are nonincreasing in x. The summands defining K's off-diagonal entries are therefore nonnegative. The theorem applies whenever the slow eigenvalue is simple.

All triplet predicates also fail directly: A<=B and t<1 imply G<0; A<=B and tu<u imply A(1+tu/2)-B(1+u/2)<0, and subtracting C(1-t) preserves strict negativity. This family can lie in the upper k branches, so it differs from the preceding corollary.

The source discusses decreasing gap laws as a class without Mpemba behavior. The contribution here is an explicit finite-temperature, all-n proof through adjacent differences. Novelty relative to all stochastic-monotonicity literature is not asserted.

## 6. Exact obstruction to the full monotone-eigenvector route

The premise of the full conjecture does **not** force an energy-monotone slow eigenvector. The following exact example disproves that intermediate implication, without claiming a Mpemba counterexample.

Take beta_b=1, q=(1,3/4,1/2,1/4), E_i=-log q_i and pi_b=q/(5/2). The downward rates in order (12,13,14,23,24,34) are

(1/4, 1/16, 3/4, 5/4, 6, 9).

The full column generator is

M = [ -13/32,  1/4,    1/16,  3/4;
        3/16, -37/12,  5/4,   6;
        1/32,  5/6,  -93/16,  9;
        3/16,  2,     9/2,  -63/4 ].

All off-diagonal entries are strictly positive; columns sum to zero; M_ij q_j=M_ji q_i exactly. The four triplets have:

| Triplet | H | G | F |
|---|---:|---:|---:|
| 123 | 7/8 | -77/64 | -317/192 |
| 124 | 101/16 | -87/16 | -199/32 |
| 134 | 309/32 | -69/8 | -167/16 |
| 234 | 155/12 | -5 | -14/3 |

Thus all four predicates fail strictly, all in the upper branch.

The characteristic polynomial is

x(1536x^3+38480x^2+175584x+100755)/1536.

The cubic has one root each in (-20,-19), (-6,-5), and (-67/100,-66/100): substitution gives opposite signs on each pair of endpoints, and three disjoint sign-changing intervals exhaust the three roots. Therefore its simple slow root x lies in (-.67,-.66).

The first column of adj(x I-M^T) is a left-eigenfunction w at that root, with entries

w1=x^3+1183x^2/48+10031x/96+6717/256;
w2=x^2/4+1333x/192+6717/256;
w3=x^2/16+467x/96+6717/256;
w4=3x^2/4+559x/64+6717/256.

Its adjacent differences simplify exactly to

w2-w1=-x(192x^2+4684x+18729)/192 >0;
w3-w2=-x(12x+133)/64 >0;
w4-w3= x(132x+743)/192 <0.

Each sign holds throughout the rational isolating interval. Hence no sign choice makes w monotone. Approximate eigenvalues, for orientation only, are 0,-.6694086130,-5.0752865526,-19.3073881677.

This blocks an attempted proof that infers a monotone slow eigenfunction from all four triplet failures. It does not block the stronger covariance-sign approach: a nonmonotone function can still have a covariance of fixed sign over the Gibbs curve.

## 7. Continuation

The useful remaining question is whether the feasible cone forces Cov_{pi(beta)}(w2,E) to have the needed sign even when w2 is nonmonotone. The upper entries of K are always nonnegative under the triplet failures, because G<=0 orders upward rates. Only downward-cumulative entries can violate the Metzler criterion. A potentially useful intermediate cone is therefore larger than the coordinatewise-positive adjacent-difference cone and depends on bath weights and the Gibbs tangent curve. No invariant such larger cone is established here.

## 8. A concise corrected four-state implication

Under the original four-triplet failure premise, all upper off-diagonal entries of K are already nonnegative. Indeed K13>=0 follows from G124<=0, K23>=0 from G234<=0, and K12 is the sum of two nonnegative terms provided by G123<=0 and G124<=0 after detailed balance.

Therefore the following three additional linear conditions are sufficient to repair the proposed implication:

1. a12>=a13;
2. a13>=a14;
3. a13+a23>=a14+a24.

They are precisely K21>=0, K31>=0, and K32>=0. Together with the original premise they imply no asymptotic cooling or heating Mpemba effect in the simple-slow-mode domain. They are sufficient, not asserted necessary or optimally weak. This corrected sufficient theorem has no energy-spacing restriction and requires only the original rates.

The verifier agent independently checked sections 1–5 and confirmed the algebra and scope. This feedback is secondary to the written proofs; it is not a substitute for them.
