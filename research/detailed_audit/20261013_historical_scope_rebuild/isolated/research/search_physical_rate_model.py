"""Reproducible, bounded numerical discovery; does not certify absence."""
import json
import time
import argparse
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from physical_rate_model import Evolution

ROOT = Path(__file__).resolve().parent
seed = 12102026
evaluations = 0
best = 1e10
history = []
start = time.monotonic()
profile = "nominal"


def perturb(x, factors):
    y = x * factors
    y[-2:] = 1 + (x[-2:]-1)*factors[-2:]
    return y


def box_rate_bound(x):
    g1,g2,g3,hA,hB,hX,_,_ = x
    maximum_activation = max(hA+g1+g2+g3, hB+g2,
                             hX+g1, hX+g1+g2, hX+g2+g3, hX+g3)
    return np.exp(1.05*maximum_activation - .95*min(hA,hB,hX))


def score(x):
    evo = Evolution(x)
    initial = evo.distances(0)[:, 0]
    if min(initial) <= .011:
        return 10 + .011 - min(initial)
    deltas = [.001, .0031622776601683794, .01]
    times = [evo.hitting(d, 0) for d in deltas]
    return float(max(evo.distances(.95*np.array(times))[1] / deltas - 1))


def objective(x):
    global evaluations, best
    evaluations += 1
    if sum(x[:3]) > 4.76 or x[7] < x[6] + .3:
        return 10 + max(sum(x[:3])-4.76, x[6]+.3-x[7])
    evo = Evolution(x)
    rates = evo.M[~np.eye(4, dtype=bool)]
    ratio = max(rates)/min(rates)
    limit_ratio = ratio if profile == "nominal" else box_rate_bound(x)
    limit = 850 if profile == "nominal" else 995
    if limit_ratio > limit:
        return 10 + np.log(limit_ratio/limit)
    # Sign of this objective tests a 5% time saving on three thresholds.
    value = score(x)
    if profile == "robust_screen":
        # Nominal and all single-coordinate extremes: discovery only.
        for j in range(8):
            for sign in [-1,1]:
                f=np.ones(8);f[j] += .05*sign
                value=max(value,score(perturb(x,f)))
    if value < best:
        best = float(value)
        history.append(dict(evaluation=evaluations, value=best, parameters=x.tolist()))
    return value


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('--profile', choices=['nominal','robust_screen'], default='nominal')
    parser.add_argument('--iterations', type=int, default=200)
    args=parser.parse_args()
    profile=args.profile
    # Interior search leaves room for a subsequent 5% neighborhood.
    bounds = [(.11, 2), (.11, 2), (.11, 4.5),
              (2.11, 5.7), (2.11, 5.7), (2.11, 5.7),
              (1.4, 5), (2, 9.5)]
    result = differential_evolution(objective, bounds, seed=seed,
                                    maxiter=args.iterations, popsize=12, polish=True,
                                    updating="immediate", tol=1e-7)
    receipt = dict(model="two-basin shared barrier excess", seed=seed, profile=profile,
                   method="scipy differential_evolution", bounds=bounds,
                   iterations=int(result.nit), evaluations=evaluations,
                   elapsed_seconds=time.monotonic()-start,
                   score=float(result.fun), history=history,
                   screening_only=True, result=Evolution(result.x).metrics())
    suffix='' if profile=='nominal' else '_'+profile
    (ROOT/('physical_rate_search'+suffix+'.json')).write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['history','result']}, indent=2))
    print(json.dumps({k:v for k,v in receipt['result'].items() if k not in ['thresholds','hitting_times_tau']}, indent=2))
