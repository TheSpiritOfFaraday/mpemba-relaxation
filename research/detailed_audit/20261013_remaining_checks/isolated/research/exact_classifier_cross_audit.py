"""Reproducible cross-audit: exact Sturm decisions vs spectral/thermal geometry.

The exact implementation is the acceptance authority. The independent numerical
calculation uses 70-digit symmetric diagonalization and the separate derivative
curve classifier. Numerical comparisons with unresolved boundary signs are
reported as unresolved, never silently counted as agreement.
"""
import json, random, time
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as s
from exact_commensurate_classifier import classify, PAIRS, require
from four_state_classifier import classify_mode

SEED = 202610112133
rng = random.Random(SEED)
mp.mp.dps = 70
start = time.time()
rows = []
outpath = Path(__file__).with_suffix('.json')

def high(v):
    v = s.Rational(v)
    return mp.mpf(int(v.p))/int(v.q)

for sample in range(100):
    base = rng.choice([s.Rational(1,2), s.Rational(2,3), s.Rational(3,4)])
    m = [0] + sorted(rng.sample(range(1,9), 3))
    q = s.Matrix([base**mi for mi in m]); q /= sum(q)
    D = s.diag(*q)
    if sample < 50:
        rates = [s.Rational(rng.randint(1,99),rng.choice([1,10,100])) for _ in PAIRS]
        family = 'independent rational rates'
    else:
        w = s.Matrix([rng.randint(-100,100) for _ in range(4)])
        w -= s.ones(4,1)*(q.T*w)[0]
        norm = (w.T*D*w)[0]
        require(norm > 0, 'nonconstant test mode')
        # Alternate well separated and deliberately close slow spectra.
        epsilon = s.Rational(1,10 if sample%2 else 10**8)/(1+max(abs(w[i]*w[j])/norm for i,j in PAIRS))
        M = q*s.ones(1,4)-s.eye(4)+epsilon*D*w*w.T/norm
        rates = [M[i,j] for i,j in PAIRS]
        family = 'prescribed rational mode'
    exact = classify(base,m,rates)
    M = s.zeros(4)
    for (i,j),rate in zip(PAIRS,rates):
        M[i,j] = rate; M[j,i] = rate*q[j]/q[i]
    for j in range(4):M[j,j] = -sum(M[i,j] for i in range(4) if i!=j)
    mq = [high(v) for v in q]
    symmetric = mp.matrix([[high(M[i,j])*mp.sqrt(mq[j]/mq[i]) for j in range(4)] for i in range(4)])
    eigenvalues,V = mp.eigsy(symmetric)
    mode = [V[i,2]/mp.sqrt(mq[i]) for i in range(4)]
    scale = max(map(abs,mode)); mode = [v/scale for v in mode]
    E = [-mi*mp.log(high(base)) for mi in m]
    floating = classify_mode([float(v) for v in E], [float(v) for v in mode], tol=1e-9)
    compared = {}
    for side in ['cooling','heating']:
        observed = floating[side]
        if isinstance(observed,(bool,np.bool_)):
            require(bool(observed)==exact['ordinary_'+side], f'sample {sample}: {side} disagreement')
            compared[side] = 'agrees'
        else:compared[side] = 'numerically unresolved'
    row = dict(sample=sample,family=family,base=str(base),exponents=m,rates=list(map(str,rates)),
        exact_cooling=exact['ordinary_cooling'],exact_heating=exact['ordinary_heating'],
        comparison=compared,slow_gap=str(eigenvalues[2]-eigenvalues[1]),
        algebraic_sign_refinements=exact['algebraic_sign_refinements'])
    rows.append(row)
    if (sample+1)%10==0:print('completed',sample+1,'elapsed',round(time.time()-start,2),flush=True)
    out = dict(status='running',seed=SEED,precision_digits=70,samples_completed=len(rows),rows=rows)
    outpath.write_text(json.dumps(out,indent=2)+'\n')
out.update(status='finished',verified=True,elapsed_seconds=time.time()-start,
    disagreements=0,unresolved_comparisons=sum(v=='numerically unresolved' for r in rows for v in r['comparison'].values()),
    cooling_cases=sum(r['exact_cooling'] for r in rows),heating_cases=sum(r['exact_heating'] for r in rows))
outpath.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
