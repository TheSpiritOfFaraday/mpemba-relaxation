"""Additional exact bounds needed by the higher-dimensional cloning proof.

Reads the rational generator formulas already proved in
quarter_limit_certificate.json. This is an explicit dependency on that
continuous-family certificate, not an independent replacement for it.
No floating-point comparison is used for acceptance.
"""
from pathlib import Path
from math import comb
import json
import sympy as s

root = Path(__file__).parent
base = json.loads((root / 'quarter_limit_certificate.json').read_text())
e, t = s.symbols('e t', positive=True)

def require(ok, message):
    if not ok:
        raise ValueError(message)

def bernstein(poly):
    p = s.Poly(s.expand(poly.subs(e, t / 30)), t)
    n = p.degree()
    return [s.factor(sum(p.nth(j) * s.Rational(comb(k, j), comb(n, j))
                        for j in range(k + 1))) for k in range(n + 1)]

receipts = []
for i in range(3):
    rate = s.sympify(base['generator'][i][3], locals={'e': e})
    num, den = s.fraction(s.cancel(2 - rate))
    nc, dc = bernstein(num), bernstein(den)
    if all(v < 0 for v in dc):
        nc, dc = [-v for v in nc], [-v for v in dc]
    require(all(v > 0 for v in nc + dc), 'a_i4 < 2 on full interval')
    receipts.append(dict(rate=f'a_{i+1}4', formula=str(rate),
                         numerator_bernstein=list(map(str, nc)),
                         denominator_bernstein=list(map(str, dc))))

# Independent symbolic checks of the three low-pair/clone H expressions.
r = s.Symbol('r', positive=True)
Q = s.Rational(13, 5) - e/3 + e**2
q = [1/Q, (1-e/3)/Q, s.Rational(3, 5)/Q, e**2*r/Q]
kcross = [1, 2, 3]
low = {(0, 1): 100, (0, 2): 10, (1, 2): 11}
expected = [(-1-98*e/3-e**2*r/2)/Q,
            (-s.Rational(24, 5)-e**2*r)/Q,
            (-s.Rational(21, 5)+3*e-e**2*r/2)/Q]
for ((i, j), kij), target in zip(low.items(), expected):
    A, B, C = q[i]*kcross[i], q[j]*kcross[j], q[i]*kij
    F = A*(1+q[3]/(2*q[i]))-B*(1+q[3]/(2*q[j]))-C*(1-q[j]/q[i])
    require(s.factor(F-target) == 0, 'low-pair clone F identity')

require(all(v > 0 for v in bernstein(3-Q)), 'Q < 3')
require(s.Rational(3, 5) < 1-s.Rational(1, 90), 'minimum low raw weight')
require(-s.Rational(21, 5)+s.Rational(3, 30) < -1, 'third F numerator <= -1')
require(s.Rational(24, 5) > 1, 'second F numerator <= -1')
require(6*s.Rational(1, 3)-s.Rational(5, 3) == s.Rational(1, 3),
        'six-epsilon-squared repair margin')

# Explicit Davis--Kahan/resolvent constants for n=5, r=(2/3,1/3).
perturb_coefficient = 6*200
epsilon_cap = s.Rational(1, 200)
perturb_cap = perturb_coefficient*epsilon_cap**2
gap = s.Rational(1, 3)-2*perturb_cap
require(gap == s.Rational(41, 150) and gap > s.Rational(1, 4), 'n5 gap')
weighted_coefficient = s.Rational(3, 2)*perturb_coefficient/(s.Rational(1, 3)-perturb_cap)
require(weighted_coefficient < 6000, 'weighted eigenvector error')
require(6000*3 == 18000, 'coordinate error using qmin >= epsilon^2/9')

out = dict(verified=True, interval='0<epsilon<=1/30',
           input_dependency='quarter_limit_certificate.json exact generator formulas',
           rate_upper_bound='2', rate_bound_certificates=receipts,
           H_low_pair_clone_F=list(map(str, expected)),
           H_low_pair_clone_F_upper_bound='-1/3',
           split_F_deficit_upper_bound='(5/3)*epsilon^2',
           perturbation='6*epsilon^2*H', repaired_F_upper_bound='-(1/3)*epsilon^2',
           n5_clone_fractions=['2/3', '1/3'], n5_epsilon_cap=str(epsilon_cap),
           H_norm_upper_bound='200', n5_slow_gap_lower_bound=str(gap),
           n5_weighted_error_upper_bound='6000*epsilon^2',
           n5_coordinate_overlap_error_upper_bound='18000*epsilon',
           scope='Exact extra rate bounds and algebra/constants; analytic cloning, perturbation and thermal-limit arguments are in the companion audit.')
(root / 'clone_rate_bound_certificate.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rate_bound_certificates'}, indent=2))
