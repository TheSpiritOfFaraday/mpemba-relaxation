# Adversarial search and exact certificate

## Executed searches

`python3 research/adversarial_search.py` ran 20,000,000 draws with seed 934256. Adjacent energy gaps were independent `10**uniform(-3,1)`. Six downward rates were independent `10**uniform(-5,5)` and rescaled by their maximum; upward rates followed detailed balance at bath inverse temperature 1. Of those draws, 2,020,118 passed the four original branch predicates. No accepted numerical witness passed the script's sign-margin and spectral-gap filters on the 18-point initial-inverse-temperature grid `{0,.0001,.01,.1,.25,.5,.75,.9,.99,1.01,1.1,1.5,2,3,5,10,100,1000}`. This is a statement about this finite numerical search, not a no-effect theorem.

`python3 research/adversarial_signed.py` used seed 742908 and sampled signed reversible perturbations. A bath-perpendicular Householder projection selected the largest eigenvalue on the three-dimensional nonstationary subspace; a subsequent positive reset component can restore positive off-diagonal rates without changing its eigenvectors. Alternating sampling families multiplied signed rates by the lower state's bath weight. Energy gaps were `10**uniform(-3,.7)`; signed scales were normal variables times `10**uniform(-3,3)`. The first candidate was reached after 9,950,000 draws, with 573,427 accepted by all eight linear F/G inequalities. It had a cooling anomaly at beta=0 with a strict derivative margin, allowing positive beta witnesses by continuity. Data are in `adversarial_signed_candidate.json`; the executed search was stopped at that first witness.

## Rational simplification

The candidate simplified to unnormalized bath weights

`q=(1,199/200,2/5,1/400)`

and signed downward perturbation rates, in pair order 12,13,14,23,24,34,

`(-24,3,-1/10,200,1/25,1/25)`.

Adding reset rates `25*q_i` for each destination i gives strictly positive downward rates

`(1,28,249/10,1799/8,4983/200,251/25)`.

Both F and G of every triplet are strictly negative, so the original branch conditions fail. The initial inverse temperatures `1/1000` and `1/100` give a strict cooling witness. With the left eigenvector normalized to `w_1=1`, the certified margin `|b_c|-|b_h|` exceeds `0.00010340762800011847331527`.

## Exact certificate

Run `python3 research/adversarial_certificate.py`. It uses SymPy rational/integer arithmetic to verify positive off-diagonal entries, column sums, detailed balance, the original branch predicates, the cubic characteristic polynomial, three disjoint negative simple root intervals, cofactor eigenvector identities modulo the cubic, zero bath expectation, and a positive normalization denominator. Interval operations are elementary exact rational operations. Bounds for `q_i**(1/1000)` and `q_i**(1/100)` are rational and certified by integer powers. Decimal approximations only choose prospective bounding rationals; every accepted bound and final inequality is checked exactly.

The right eigenvector is `v=diag(q)w/(sum_i q_i*w_i**2)`, so `w^T v=1`. The eigenvalue is the largest of the three negative roots. The JSON contains the exact bounds; the log contains readable decimal approximations.

`adversarial_rational_mode_certificate.py` certifies an auxiliary variant with exactly rational w and lambda, at the price of more complicated rates. It retains the same q, initial temperatures, and strict predicate failure. Its w and lambda are

`w=(1,-759/1000,-48957/80000,-1/250)`, `lambda=-4966671/250000`.

The companion is useful corroboration and is not needed for the primary certificate.

## Mechanism

The very rapid 2–3 transitions make those two left-mode components close, while state 1 has the opposite sign. The fourth state's slow-mode component is close to zero. At very high initial temperatures the Gibbs curve changes the contribution of state 4 strongly enough that the slow-mode overlap initially becomes more negative as inverse temperature increases. It subsequently rises to zero at the bath temperature. Thus a hotter finite-temperature initial state has a smaller magnitude of the slow-mode amplitude. The original triplet inequalities constrain local rates but miss the global slow-mode orientation responsible for this turning behavior.

## Exact family

For every real `alpha>24`, the downward rates

`(alpha-24, alpha+3, alpha-1/10, (199/200)*alpha+200, (199/200)*alpha+1/25, (2/5)*alpha+1/25)`

give a strictly positive counterexample with the same q, slow eigenvectors, and finite-temperature witness. They satisfy

`M(alpha)=M(25)+(alpha-25)*(q*1^T-(959/400)*I)`.

On nonstationary right eigenvectors, the reset term acts as `-(959/400)*I`; consequently all three nonstationary eigenvalues shift equally, and simplicity/order are preserved. Each F and G is invariant under this reset addition (verified algebraically in the primary certificate). Positivity follows directly from alpha>24. Hence the example is an exact one-parameter family. Strict predicate margins, spectral separation, and witness margin also imply an open neighborhood of counterexamples by continuity.
