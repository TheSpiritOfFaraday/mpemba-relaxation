"""Audit manually transcribed sharpness formulas and printed enclosures.

This supplements the analytic review and prior full certificate replay. It
does not turn receipt data into an independent proof of the physical family.
"""
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
checks = []


def require(condition, name):
    if not condition:
        raise ValueError(name)
    checks.append(name)


def zero(value):
    return s.cancel(value) == 0


e = s.Symbol('e', positive=True)
R = s.Rational
raw = s.Matrix([1, 1-e/3, R(3, 5), e**2])
w = s.Matrix([1, -R(4, 5), -R(1, 3)-4*e/9+5*e**3/3, -e])
P = 15*e**4-5*e**3+69*e**2-20*e+68
T = 15*e**3-4*e-12
S = 75*e**3-20*e+21
rates = [s.Integer(1),
         (150*e**6+640*e**5+190*e**4+459*e**3+138*e**2-252*e-1224)/(T*P),
         (30*e**3+125*e**2-42*e+240)/(3*P),
         5*(3-e)*(225*e**5-185*e**4+425*e**3-113*e**2-60*e+204)/(S*P),
         5*(3-e)*(9*e**2-2*e+16)/(3*P), 3*(9*e**2-2*e+16)/P]
M = s.zeros(4)
for (i, j), rate in zip(itertools.combinations(range(4), 2), rates):
    M[i, j], M[j, i] = rate, rate*raw[j]/raw[i]
for j in range(4):
    M[j, j] = -sum(M[i, j] for i in range(4) if i != j)
require(zero((raw.T*w)[0]), 'formal family bath orthogonality')
require(all(zero(v) for v in M.T*w+3*w), 'formal family four eigenvector equations')
require(all(zero(v) for v in s.ones(1, 4)*M), 'formal family column conservation')
require(all(zero(v) for v in M*s.diag(*raw)-s.diag(*raw)*M.T), 'formal family detailed balance')
require([v.subs(e, 0) for v in rates] == [1, R(3, 2), R(20, 17), R(15, 7), R(20, 17), R(12, 17)], 'six printed rate limits')
triplets = {}
for i, j, k in itertools.combinations(range(4), 3):
    a, b, c = M[i, k], M[j, k], M[i, j]
    t, u = raw[j]/raw[i], raw[k]/raw[j]
    triplets[(i, j, k)] = (s.cancel(a*(1+t*u/2)-b*(1+u/2)-c*(1-t)), s.cancel(a*t-b))
require(zero(triplets[(0, 1, 3)][0]), 'formal family F124 identically zero')
require(zero(triplets[(1, 2, 3)][1]), 'formal family G234 identically zero')
factor = 2*(e-3)*(15*e**2-5*e-6)/(9*P)
require(zero(rates[4]-rates[2]*raw[1]-e*factor), 'vanishing G124 factor extracted exactly')
a, b, c, d, f, g = rates
t, u, v = raw[1:]
F123 = b*t*(2+u)-d*(2*t+u)-2*a*t*(1-t)
F134 = c*u*(2+v)-g*(2*u+v)-2*b*u*(1-u)
F234 = f*u*(2*t+v)-g*t*(2*u+v)-2*d*u*(t-u)
require(all(zero(z) for z in [F123-2*t*triplets[(0, 1, 2)][0],
                              F134-2*u*triplets[(0, 2, 3)][0],
                              F234-2*t*u*triplets[(1, 2, 3)][0]]), 'cleared F expressions preserve the original signs')
continuous_checks = [P, -T, S, *rates[1:], d-b*t, -F123, -F134, -F234, factor]
values = [s.cancel(x).subs(e, 0) for x in continuous_checks]
require(len(values) == 13 and all(x > 0 for x in values), 'all thirteen continuity checks positive at zero')
require(all(s.denom(s.cancel(x)).subs(e, 0) != 0 for x in continuous_checks), 'all thirteen rational denominators nonzero at zero')

# Compare the manuscript's explicit U with the independently generated basis.
q = raw/sum(raw)
D = s.diag(*q)
r = w[2]
U = s.Matrix([[-R(3, 5)-(r-1)/3, e**2*(5*e-4)/9],
              [(r-1)/(3*(1-e/3)), -5*e**2*(e+1)/(9*(1-e/3))], [1, 0], [0, 1]])
require(all(zero(x) for x in q.T*U), 'printed fast basis has zero bath mean')
require(all(zero(x) for x in w.T*D*U), 'printed fast basis is orthogonal to slow mode')
cert = json.loads((ROOT/'research/quarter_limit_certificate.json').read_text())
computed_U = s.Matrix([[s.sympify(x, locals={'e': e}) for x in row] for row in cert['U']])
require(all(zero(x) for x in computed_U-U), 'printed and certificate fast bases agree exactly')
require(cert['certificate_count'] == len(cert['certificates']) == 30, 'thirty continuous sign certificates recorded')
require(all(R(x) > 0 for c in cert['certificates'].values() for key in ['numerator_coefficients', 'denominator_coefficients'] for x in c[key]), 'all recorded Bernstein coefficients strictly positive')

# Worst case of the analytic perturbation estimates: epsilon <= 1/30.
rho = R(100, 30**4)
require(R(1, 3)-2*rho > R(1, 4), 'perturbed spectral gap remains above one quarter')
require(2*100**2/(R(1, 3)-rho)**2 < 500**2, 'weighted eigenvector difference bound below 500 epsilon^4')
require(3*500**2 < 1000**2, 'coordinate difference bound below 1000 epsilon^3')
thermal = json.loads((ROOT/'research/quarter_limit_thermal_certificate.json').read_text())
lo, hi = map(R, thermal['advantage_enclosure'])
require(R('22.6171597478') < 100*lo <= 100*hi < R('22.6171597506'), 'printed finite-witness decimal endpoints are outward bounds')
require(R(thermal['strict_advantage_lower_bound']) > R(22, 100), 'strict finite-witness lower bound exceeds 22 percent')

files = ['formal/Mpemba/SharpFamily.lean', 'formal/Mpemba/PositiveMode.lean',
         'formal/Mpemba/SharpThermal.lean', 'formal/Mpemba/SharpQuarter.lean',
         'research/quarter_limit_certificate.py', 'research/quarter_limit_certificate.json',
         'research/quarter_limit_thermal_certificate.py', 'research/quarter_limit_thermal_certificate.json',
         'research/standalone_primary_check.py', 'research/quarter_limit_family.tex',
         'research/finite_span_amplitude.tex', 'focused_paper.tex', 'focused_appendices.tex']
out = {'passed': True, 'check_count': len(checks), 'checks': checks,
       'formal_family_positive_values_at_zero': list(map(str, values)),
       'source_hashes': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
       'scope': 'Exact formula correspondence, perturbation constants and printed bounds. Analytic positivity, spectral, continuity and thermal-limit reasoning are reviewed separately. Receipt coefficient inspection does not replace the prior full certificate replay.'}
(ROOT/'research/detailed_audit/sharpness_correspondence.json').write_text(json.dumps(out, indent=2)+'\n')
print('PASS:', len(checks), 'sharpness formula and enclosure comparisons')
