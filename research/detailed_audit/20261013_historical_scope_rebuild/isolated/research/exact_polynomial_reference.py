"""Independent rational-polynomial sign-change reference for audit scripts.

Uses SymPy's rational root counter directly, not OrderedField or numerical
isolating-box locations. Endpoints and even-multiplicity roots do not count.
"""
import sympy as s


def sign_changes_open(poly, left, right):
    if len(poly.gens) != 1 or poly.domain != s.QQ:
        raise ValueError('univariate rational polynomial required')
    if not left < right or poly.is_zero:
        raise ValueError('ordered interval and nonzero polynomial required')
    x = poly.gens[0]
    odd = s.Poly(1, x, domain=s.QQ)
    for factor, multiplicity in poly.sqf_list()[1]:
        if multiplicity % 2:
            odd *= factor
    for endpoint in (left, right):
        while odd.degree() > 0 and odd.eval(endpoint) == 0:
            odd, remainder = odd.div(s.Poly(x-endpoint, x, domain=s.QQ))
            if not remainder.is_zero:
                raise ValueError('exact endpoint factor did not divide')
    return int(odd.count_roots(left, right)) if odd.degree() > 0 else 0
