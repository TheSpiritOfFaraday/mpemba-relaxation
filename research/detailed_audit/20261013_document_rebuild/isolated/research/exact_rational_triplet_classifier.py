"""Exact full-triplet four-state classification without commensurate energies.

Input: four strictly decreasing positive rational bath weights and six positive
rational downward rates. Energies are -log(weights), bath inverse temperature 1.
Requires every original triplet predicate to fail. Uses the scalar endpoint
quadratic, exact algebraic coefficients, and Baker's logarithm independence.
No floating-point value is used to accept an input or decide a sign.

This is a mathematical algorithm with a symbolic/interval implementation, not
a Lean proof. General unconstrained systems and finite-time KL/TV are outside
its scope. See rational_log_endpoint.tex for the termination argument.
"""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from exact_commensurate_classifier import OrderedField, require

X = s.Symbol('lambda')
PAIRS = list(itertools.combinations(range(4), 2))


def rational(value):
    require(not isinstance(value, float), 'use rational strings or integers, not floats')
    result = s.sympify(value)
    require(bool(result.is_Rational), 'exact rational input required')
    return s.Rational(result)


def interval_product(left, right):
    values = [a * b for a in left for b in right]
    return min(values), max(values)


def algebraic_bounds(field, value):
    if field.rational:
        v = s.Rational(value.numerator, value.denominator)
        return v, v
    interval = field.interval
    lo = s.Rational(interval.a.numerator, interval.a.denominator)
    hi = s.Rational(interval.b.numerator, interval.b.denominator)
    result = (s.Rational(0), s.Rational(0))
    for c in value.to_list():
        c = s.Rational(c.numerator, c.denominator)
        a, b = interval_product(result, (lo, hi))
        result = (a + c, b + c)
    return result


def log_ratio_bounds(ratio, terms):
    """Exact atanh-series enclosure, restricted to 1 <= ratio <= 2."""
    require(1 <= ratio <= 2 and terms > 0, 'logarithm range and series length')
    z = (ratio - 1) / (ratio + 1)
    total = s.Rational(0)
    power = z
    for k in range(terms):
        total += 2 * power / (2 * k + 1)
        power *= z * z
    remainder = 2 * power / ((2 * terms + 1) * (1 - z * z))
    return total, total + remainder


def log_prime_bounds(prime, terms):
    exponent = int(prime).bit_length() - 1
    ratio = s.Rational(prime, 2 ** exponent)
    a, b = log_ratio_bounds(ratio, terms)
    c, d = log_ratio_bounds(s.Rational(2), terms)
    return a + exponent * c, b + exponent * d


def prime_exponents(weight):
    result = {}
    for p, count in s.factorint(weight.p).items():
        result[int(p)] = int(count)
    for p, count in s.factorint(weight.q).items():
        result[int(p)] = result.get(int(p), 0) - int(count)
    return result


def classify(weights, rates, max_rounds=None):
    weights = list(map(rational, weights))
    rates = list(map(rational, rates))
    require(len(weights) == 4 and all(weights[i] > weights[i + 1] > 0 for i in range(3)),
            'four strictly decreasing positive weights')
    require(len(rates) == 6 and all(v > 0 for v in rates), 'six strictly positive rates')
    if max_rounds is not None:
        require(isinstance(max_rounds, int) and max_rounds > 0, 'positive round limit')
    q = s.Matrix(weights) / sum(weights)
    D = s.diag(*q)
    M = s.zeros(4)
    for (i, j), rate in zip(PAIRS, rates):
        M[i, j] = rate
        M[j, i] = rate * q[j] / q[i]
    for j in range(4):
        M[j, j] = -sum(M[i, j] for i in range(4) if i != j)
    require(M * D == D * M.T, 'detailed balance')
    require(s.ones(1, 4) * M == s.zeros(1, 4), 'column conservation')
    triplets = {}
    for i, j, k in itertools.combinations(range(4), 3):
        A, B, C = M[i, k], M[j, k], M[i, j]
        F = A * (1 + q[k] / (2 * q[i])) - B * (1 + q[k] / (2 * q[j])) - C * (1 - q[j] / q[i])
        G = A * q[j] / q[i] - B
        require(F <= 0 and G <= 0, 'input outside the full triplet-failure premise')
        triplets[f'{i+1}{j+1}{k+1}'] = {'F': str(F), 'G': str(G)}
    char = s.Poly(M.charpoly(X).as_expr() / X, X)
    roots = char.all_roots(radicals=False)
    require(all(root.is_real for root in roots), 'real reversible spectrum')
    alpha = roots[-1]
    multiplicity = sum(root == alpha for root in roots)
    field = OrderedField(alpha)
    require(field.sign(field.convert(alpha)) < 0, 'strictly negative slow rate')
    result = dict(
        verified=True, weights=list(map(str, weights)), downward_rates=list(map(str, rates)),
        triplets=triplets, nonzero_characteristic_polynomial=str(char.as_expr()),
        slow_eigenvalue=str(alpha), slow_multiplicity=multiplicity,
        ordinary_heating=False, strong_cooling=False, strong_heating=False,
        scope='Full-triplet premise only; asymptotic cooling/heating in KL and chi-squared, any fixed norm for simple slow mode. No general finite-time classification.',
    )
    if multiplicity > 1:
        result.update(ordinary_cooling=False, method='exact repeated-slow classification',
                      all_time_convex_divergence_reversal=False)
        return result
    tau = s.trace(M)
    sigma = (tau ** 2 - s.trace(M * M)) / 2
    a = s.ones(4, 1) / 4 - q
    h = D.inv() * a
    factors = list(map(prime_exponents, weights))
    primes = sorted(set().union(*(f.keys() for f in factors)))
    coefficients = {}
    polynomial_coefficients = {}
    for prime in primes:
        exponents = s.Matrix([f.get(prime, 0) for f in factors])
        d = (exponents - s.ones(4, 1) * sum(exponents) / 4) / 4
        A0, A1, A2 = [(h.T * M ** j * d)[0] for j in range(3)]
        polynomial = A0 * X ** 2 + (A1 - tau * A0) * X + A2 - tau * A1 + sigma * A0
        coefficient = field.convert(polynomial.subs(X, alpha))
        polynomial_coefficients[str(prime)] = str(s.expand(polynomial))
        if coefficient:
            coefficients[prime] = coefficient
    result['prime_log_coefficient_polynomials'] = polynomial_coefficients
    result['nonzero_prime_log_coefficients'] = {
        str(prime): str(field.K.to_sympy(value)) for prime, value in coefficients.items()
    }
    if not coefficients:
        result.update(ordinary_cooling=False, endpoint_sign=0,
                      method='exact algebraic cancellation of every prime-log coefficient')
        return result
    # Baker's theorem guarantees this nonzero algebraic linear combination
    # of logarithms of distinct positive primes cannot vanish.
    round_index = 0
    while max_rounds is None or round_index < max_rounds:
        bits = 16 * 2 ** round_index
        if not field.rational:
            while s.Rational(field.interval.dx.numerator, field.interval.dx.denominator) > s.Rational(1, 2 ** bits):
                field.interval = field.interval.refine()
                field.refinements += 1
        lower = upper = s.Rational(0)
        for prime, coefficient in coefficients.items():
            lo, hi = interval_product(algebraic_bounds(field, coefficient), log_prime_bounds(prime, bits))
            lower += lo
            upper += hi
        round_index += 1
        if lower > 0 or upper < 0:
            sign = 1 if lower > 0 else -1
            result.update(ordinary_cooling=sign > 0, endpoint_sign=sign,
                          endpoint_interval=[str(lower), str(upper)],
                          method='Baker prime-log independence and exact rational interval signs',
                          interval_rounds=round_index, algebraic_refinements=field.refinements)
            return result
    result.update(verified=False, ordinary_cooling=None, endpoint_sign=None,
                  method='round limit reached; unresolved, not absence', interval_rounds=round_index)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='JSON containing bath_weights and downward_rates')
    parser.add_argument('--output')
    parser.add_argument('--max-rounds', type=int)
    args = parser.parse_args()
    data = json.loads(Path(args.input).read_text())
    result = classify(data['bath_weights'], data['downward_rates'], args.max_rounds)
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end='')
