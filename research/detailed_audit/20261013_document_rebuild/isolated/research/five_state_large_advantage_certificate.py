"""Exact thermal witness for the five-state cloned, repaired family.

The generator and triplet signs are proved in high_dimension_amplitude.tex.
The n=5 perturbation estimate is independently recorded in the cloning audit:
weighted error <=6000 e^2 and hence component error <=18000 e for e<=1/200.
All thermal enclosures below use integer powers and rational arithmetic.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import time
from standalone_primary_check import approximate_overlap, require

started = time.time()
epsilon = F(1, 10**50)
q = (F(1), 1-epsilon/3, F(3, 5), 2*epsilon**2/3, epsilon**2/3)
w = (F(1), F(-4, 5), -F(1, 3)-4*epsilon/9+5*epsilon**3/3, -epsilon, -epsilon)
require(sum(a*b for a, b in zip(q, w)) == 0, 'zero bath overlap')
hot = approximate_overlap(q, w, 10000)
cold = approximate_overlap(q, w, 50)
error = 18000 * epsilon
hot_repaired = (hot[0]-error, hot[1]+error)
cold_repaired = (cold[0]-error, cold[1]+error)
require(hot_repaired[1] < 0 and cold_repaired[1] < 0, 'strict overlap signs')
lower_gain = 1 - (-hot_repaired[0])/(-cold_repaired[1])
upper_gain = 1 - (-hot_repaired[1])/(-cold_repaired[0])
require(lower_gain > F(35, 100), 'strict repaired amplitude advantage exceeds35percent')
outward_percent = ('37.6778168708', '37.6778168731')
require(F(outward_percent[0]) < 100*lower_gain and 100*upper_gain < F(outward_percent[1]),
        'outward decimal display bounds')
result = dict(
    verified=True, epsilon=str(epsilon), bath_weights=list(map(str, q)),
    unperturbed_slow_mode=list(map(str, w)), beta_hot='1/10000', beta_cold='1/50',
    beta_bath='1', perturbation='6 epsilon^2 H as in high_dimension_amplitude.tex',
    overlap_error_bound=str(error), hot_overlap=list(map(str, hot_repaired)),
    cold_overlap=list(map(str, cold_repaired)),
    repaired_gain_enclosure=list(map(str, (lower_gain, upper_gain))),
    displayed_gain_percent=[100*float(v) for v in (lower_gain, upper_gain)],
    displayed_gain_percent_are_approximations=True,
    outward_decimal_gain_percent=list(outward_percent),
    elapsed_seconds=time.time()-started,
    scope='Exact thermal enclosures plus explicit analytic spectral-perturbation bound; positivity/triplets/gap depend on the companion proof and symbolic certificate. Not a Lean formalization.',
)
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: result[k] for k in ['verified', 'displayed_gain_percent', 'elapsed_seconds', 'scope']}, indent=2))
