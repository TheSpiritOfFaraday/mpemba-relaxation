"""Exact obstruction to coordinatewise KL monotonicity on the contrast cone.

This supplies no Markov trajectory or thermal-preparation pair. It cannot
establish a finite-time heating Mpemba effect.
"""
from fractions import Fraction as Q
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


q = [Q(4, 10), Q(3, 10), Q(2, 10), Q(1, 10)]
x = [Q(3, 4), Q(1, 100), Q(1, 100)]
H = [[Q(0) for _ in x] for _ in q]
for k in range(3):
    for i in range(k + 1):
        H[i][k] = Q(1)
    H[k + 1][k] = -sum(q[:k + 1]) / q[k + 1]
r = [1 + sum(H[i][k] * x[k] for k in range(3)) for i in range(4)]
require(r == [Q(177, 100), Q(1, 50), Q(39, 40), Q(91, 100)],
        'Unexpected likelihood ratios')
require(all(z > 0 for z in r), 'Probability boundary reached')
require(sum(a * b for a, b in zip(q, r)) == 1, 'Normalization failed')
require(sum(q[i] * H[i][1] for i in range(4)) == 0, 'Tangent not centered')
# d KL(q*r || q)/dx_2 = (1/10) log(r_1^4 r_2^3/r_3^7).
product = r[0]**4 * r[1]**3 / r[2]**7
require(0 < product < 1, 'Negative KL derivative not established')
eta2 = sum(q[i] * H[i][1]**2 for i in range(4))
chi_derivative = 2 * eta2 * x[1]
require(chi_derivative > 0, 'Positive chi-square derivative not established')
receipt = {
    'verified': True,
    'bath_weights': list(map(str, q)),
    'positive_contrast_coordinates': list(map(str, x)),
    'likelihood_ratios': list(map(str, r)),
    'exp_of_ten_times_KL_derivative': str(product),
    'KL_derivative_strictly_negative': True,
    'chi_square_derivative': str(chi_derivative),
    'scope': 'The KL function is not coordinatewise increasing on the entire '
             'positive weighted-prefix contrast cone. No common semigroup '
             'or Gibbs preparation pair is constructed; this is not a '
             'finite-time heating counterexample.'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
