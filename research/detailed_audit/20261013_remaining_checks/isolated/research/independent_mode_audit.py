"""Independent exact audit; no imports from discovery/certificate scripts.

Run from the workspace: python3 research/independent_mode_audit.py
SymPy supplies exact rational matrix arithmetic and Sturm root counts.
"""
from itertools import combinations
import sympy as s

R = s.Rational
q = s.Matrix([1, R(199, 200), R(2, 5), R(1, 400)])
w = s.Matrix([1, -R(759, 1000), -R(48957, 80000), -R(1, 250)])
slow = -R(4966671, 250000)
rates = [R(1), R(28), R(249, 10), R(330154675693, 1470375000),
         R(4983, 200), R(1527327066, 151990625)]
M = s.zeros(4)
for (i, j), a in zip(combinations(range(4), 2), rates):
    M[i, j], M[j, i] = a, a*q[j]/q[i]
for j in range(4):
    M[j, j] = -sum(M[i, j] for i in range(4) if i != j)

assert all(sum(M[:, j]) == 0 for j in range(4))
assert all(M[i, j] > 0 and M[i, j]*q[j] == M[j, i]*q[i]
           for i in range(4) for j in range(4) if i != j)
assert M.T*w == slow*w
assert (q.T*w)[0] == 0
norm = sum(q[i]*w[i]**2 for i in range(4))
v = s.diag(*q)*w/norm
assert (w.T*v)[0] == 1
assert M*v == slow*v
print('Exact left/right normalization:', norm)

for i, j, k in combinations(range(4), 3):
    A, B, C = M[i, k], M[j, k], M[i, j]
    t, u = q[j]/q[i], q[k]/q[j]
    kappa = (A+B)/(C*(1+t))
    F = A*(1+t*u/2)-B*(1+u/2)-C*(1-t)
    G = A*t-B
    predicate = ((kappa >= 1 and F > 0) or (kappa <= 1 and G > 0))
    assert not predicate
    assert F < 0 and G < 0
    print('Triplet', (i+1, j+1, k+1), 'kappa,F,G:', kappa, F, G)

x = s.Symbol('x')
actual = s.factor(M.charpoly(x).as_expr())
quadratic = (2530028648200000000*x**2 + 1002356769715761942000*x
             + 50981800562345565497767)
expected = x*(250000*x+4966671)*quadratic/632507162050000000000000
assert s.expand(actual-expected) == 0
poly = s.Poly(quadratic, x)
assert poly.count_roots(-337, -336) == 1
assert poly.count_roots(-60, -59) == 1
assert -59 < slow < 0
print('Exact characteristic polynomial:', actual)
print('Other roots: one in (-337,-336), one in (-60,-59).')
print('Simple slow eigenvalue:', slow)
print('ALL INDEPENDENT EXACT SPECTRAL AND PREDICATE CHECKS PASSED')
