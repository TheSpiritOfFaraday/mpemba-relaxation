"""Alternative exact proof of the PRIMARY four-state counterexample.

Python standard library only; no file inputs, SymPy, floating point, cofactor
polynomials, eigenvalue approximations accepted on trust, or prior certificates.
Run anywhere: python3 -I -S /absolute/path/standalone_primary_check.py

Spectral counts: exact LDL congruence of diag(q)*(L - shift*I).
Mode control: weighted residual and a certified spectral gap, not root/cofactor
interval evaluation. Thermal weights: integer binary search for nth roots.
The printed JSON is a receipt; the proof checks run before it is printed.
This is algorithmic independence, not independent external authorship or review.
"""
from fractions import Fraction as F
from itertools import combinations
import json
import sys

Q = tuple(map(F, ('1', '199/200', '2/5', '1/400')))
RATES = tuple(map(F, ('1', '28', '249/10', '1799/8', '4983/200', '251/25')))
R = F('-19.86836952145135')
W_PREFIX = tuple(map(F, ('1', '-0.75893729814493', '-0.61211930644096')))
ETA = F(1, 10**8)
GAP = F(39)


class VerificationError(Exception):
    pass


def require(condition, message):
    # Deliberately not assert: these checks also run under python -O.
    if not condition:
        raise VerificationError(message)


def row_generator(q, rates):
    require(len(q) == 4 and len(rates) == 6, 'input dimensions')
    require(all(q[i] > q[i+1] > 0 for i in range(3)), 'ordered positive weights')
    require(all(a > 0 for a in rates), 'strict rate positivity')
    L = [[F(0) for _ in q] for _ in q]
    for (i, j), downhill in zip(combinations(range(4), 2), rates):
        L[j][i] = downhill                 # from high j to low i
        L[i][j] = downhill * q[j] / q[i]   # reverse thermal rate
    for i in range(4):
        L[i][i] = -sum(L[i])
    require(all(sum(row) == 0 for row in L), 'row conservation')
    require(all(q[i]*L[i][j] == q[j]*L[j][i]
                for i in range(4) for j in range(4)), 'detailed balance')
    return L


def positive_inertia(L, q, shift):
    """Number of eigenvalues of L above shift, via exact symmetric elimination."""
    A = [[q[i]*(L[i][j] - (shift if i == j else 0))
          for j in range(4)] for i in range(4)]
    require(all(A[i][j] == A[j][i] for i in range(4) for j in range(4)),
            'symmetric pencil')
    pivots = []
    for k in range(4):
        d = A[k][k]
        require(d != 0, 'zero LDL pivot; no inference permitted')
        pivots.append(d)
        for i in range(k+1, 4):
            for j in range(k+1, 4):
                A[i][j] -= A[i][k]*A[k][j]/d
    return sum(d > 0 for d in pivots), pivots


def thermal_root(q, n, scale=10**12):
    """Enclose q**(1/n) without evaluating a real power or using a root library."""
    require(0 < q <= 1 and n > 0, 'thermal root domain')
    target = q.numerator * scale**n
    denominator = q.denominator
    lo, hi = 0, scale+1
    # Invariant: lo**n*den <= target < hi**n*den.
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**n * denominator <= target:
            lo = mid
        else:
            hi = mid
    lower, upper = F(lo, scale), F(hi, scale)
    require(lower**n <= q < upper**n, 'integer-power root enclosure')
    return lower, upper


def approximate_overlap(q, w, n):
    """A scalar enclosure of w.pi(1/n), with signs treated explicitly."""
    roots = [thermal_root(qi, n) for qi in q]
    numerator_lo = sum(wi*(lo if wi >= 0 else hi)
                       for wi, (lo, hi) in zip(w, roots))
    numerator_hi = sum(wi*(hi if wi >= 0 else lo)
                       for wi, (lo, hi) in zip(w, roots))
    denominator_lo = sum(lo for lo, hi in roots)
    denominator_hi = sum(hi for lo, hi in roots)
    require(0 < denominator_lo <= denominator_hi, 'partition denominator')
    corners = [a/b for a in (numerator_lo, numerator_hi)
               for b in (denominator_lo, denominator_hi)]
    return min(corners), max(corners)


def verify(q=Q, rates=RATES, prefix=W_PREFIX, trial=R, hot_n=1000, cold_n=100):
    L = row_generator(q, rates)
    require(hot_n > cold_n > 1, 'finite cooling temperature ordering')
    triplets = []
    for i, j, k in combinations(range(4), 3):
        # a_ik=L[k][i], a_jk=L[k][j], a_ij=L[j][i].
        A, B, C = L[k][i], L[k][j], L[j][i]
        t, u = q[j]/q[i], q[k]/q[j]
        kappa = (A+B)/(C*(1+t))
        F_value = A*(1+t*u/2)-B*(1+u/2)-C*(1-t)
        G_value = A*t-B
        # Check the published branches directly; no cone-equivalence lemma used.
        upper_pass = kappa >= 1 and F_value > 0
        lower_pass = kappa <= 1 and G_value > 0
        require(not (upper_pass or lower_pass), 'published triplet predicate')
        triplets.append({'states': [i+1, j+1, k+1], 'kappa': str(kappa),
                         'F': str(F_value), 'G': str(G_value)})

    counts = []
    for shift, expected in [(-337, 4), (-336, 3), (-60, 3),
                            (-59, 2), (-20, 2), (-19, 1)]:
        count, pivots = positive_inertia(L, q, F(shift))
        require(count == expected, 'spectral inertia count')
        counts.append({'shift': shift, 'positive_count': count,
                       'pivots': list(map(str, pivots))})
    # Exactly one root in each negative bracket; the remaining eigenvalue is 0.
    require(-20 < trial < -19, 'trial belongs to slow interval')
    require(trial - F(-59) > GAP, 'distance to the two faster modes')
    w = list(prefix) + [-sum(qi*wi for qi, wi in zip(q[:3], prefix))/q[3]]
    require(len(w) == 4 and sum(qi*wi for qi, wi in zip(q, w)) == 0,
            'trial vector has zero bath mean')
    residual = [sum(L[i][j]*w[j] for j in range(4))-trial*w[i]
                for i in range(4)]
    residual_squared = sum(qi*ri**2 for qi, ri in zip(q, residual))
    require(residual_squared < ETA**2, 'weighted eigenvector residual')
    require(sum(qi*wi**2 for qi, wi in zip(q, w)) > (ETA/GAP)**2,
            'nonzero slow projection')
    require(min(q) >= F(1, 400), 'coordinate bound from bath weights')
    # Weighted spectral theorem: z=P_slow w satisfies ||z-w||_q<ETA/GAP.
    # Then |(z-w).p| <= 20*ETA/GAP for ANY probability vector p.
    error = 20*ETA/GAP
    hot = approximate_overlap(q, w, hot_n)
    cold = approximate_overlap(q, w, cold_n)
    hot_true = hot[0]-error, hot[1]+error
    cold_true = cold[0]-error, cold[1]+error
    require(hot_true[1] < 0 and cold_true[1] < 0, 'both true overlaps negative')
    margin = hot_true[0]-cold_true[1]
    require(margin > F(1, 10000), 'strict slow-overlap advantage')
    require(w[0]-error > 0, 'positive first coordinate of true slow projection')
    normalized_margin = margin/(w[0]+error)
    require(normalized_margin > F(1, 10000), 'advantage with first coordinate one')
    # Conservatively rounded decimal brackets, still exact rationals.
    report_hot = F('-0.093745434'), F('-0.093745422')
    report_cold = F('-0.093848841'), F('-0.093848830')
    require(report_hot[0] < hot_true[0] <= hot_true[1] < report_hot[1],
            'displayed hot enclosure')
    require(report_cold[0] < cold_true[0] <= cold_true[1] < report_cold[1],
            'displayed cold enclosure')
    return {'verified': True, 'method': 'rational LDL inertia and weighted residual',
            'q': list(map(str, q)), 'down_rates': list(map(str, rates)),
            'column_generator': [[str(L[j][i]) for j in range(4)] for i in range(4)],
            'inverse_temperatures': [str(F(1, hot_n)), str(F(1, cold_n))],
            'left_trial': list(map(str, w)), 'trial_eigenvalue': str(trial),
            'triplets': triplets, 'inertia': counts,
            'residual_squared': str(residual_squared), 'residual_bound': str(ETA),
            'spectral_separation_bound': str(GAP), 'overlap_error_bound': str(error),
            'hot_projection_overlap': list(map(str, hot_true)),
            'cold_projection_overlap': list(map(str, cold_true)),
            'strict_margin_lower_bound': str(margin),
            'first_coordinate_one_margin_bound': str(normalized_margin),
            'reported_hot_enclosure': list(map(str, report_hot)),
            'reported_cold_enclosure': list(map(str, report_cold)),
            'normalization': 'z=P_slow w; v=diag(q)z/(z^T diag(q)z)',
            'scope': 'Primary asymptotic cooling counterexample; not finite crossing times'}


if __name__ == '__main__':
    try:
        print(json.dumps(verify(), indent=2))
    except VerificationError as exc:
        print('VERIFICATION FAILED: '+str(exc), file=sys.stderr)
        raise SystemExit(1)
