"""Recheck the final exploratory KL candidate by a direct matrix exponential.

This is a high-precision cross-check, not an interval certificate. Decimal
values in the saved candidate are treated as exact input decimals; equilibrium
is recomputed from the energies rather than copied from floating-point output.
"""
from pathlib import Path
import json
import mpmath as mp

ROOT = Path(__file__).resolve().parent
candidate = json.loads((ROOT / 'transient_heating_optimized_candidate.json').read_text())
mp.mp.dps = 100
E = mp.matrix([mp.mpf(str(v)) for v in candidate['E']])
q = mp.matrix([mp.exp(-v) for v in E])
q /= sum(q)
K = mp.matrix([[mp.mpf(str(v)) for v in row] for row in candidate['K']])
M = mp.matrix(4)
for j in range(4):
    for i in range(4):
        if i != j:
            M[i, j] = q[i] * K[i, j]
    M[j, j] = -sum(M[i, j] for i in range(4) if i != j)
beta = mp.mpf(str(candidate['beta']))
t = mp.mpf(str(candidate['time']))
p0 = mp.matrix([mp.exp(-beta * v) for v in E])
p0 /= sum(p0)
mean_energy = sum(p0[i] * E[i] for i in range(4))
dp = mp.matrix([-p0[i] * (E[i] - mean_energy) for i in range(4)])
T = mp.expm(t * M)
pt, dt = T * p0, T * dp
terms = [dt[i] * mp.log(pt[i] / q[i]) for i in range(4)]
derivative = sum(terms)
result = {
    'precision': mp.mp.dps,
    'source_candidate': candidate,
    'derivative': mp.nstr(derivative, 90),
    'score': mp.nstr(derivative / sum(abs(v) for v in terms), 90),
    'positive': bool(derivative > 0),
    'scope': 'Independent high-precision direct matrix exponential; not interval certification',
}
(ROOT / 'transient_heating_final_high_precision.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'source_candidate'}, indent=2))
