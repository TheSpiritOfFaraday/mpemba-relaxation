"""Search the continuous uncertainty box and threshold band for a failure.

A positive minimum is numerical evidence only, never a proof of absence.
"""
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
from physical_rate_model import Evolution
from check_physical_rate_model import perturbed

ROOT=Path(__file__).resolve().parent
source=json.loads((ROOT/'physical_rate_refined.json').read_text())
base=source['result']['parameters']


def objective(z):
    e=Evolution(perturbed(base,z[:8]))
    d=10**z[8]
    tc=e.hitting(d,0);th=e.hitting(d,1)
    return 1-th/tc


if __name__=='__main__':
    start=time.monotonic()
    runs=[]
    for seed in [12102031,12102032,12102033]:
        result=differential_evolution(objective,[(-1,1)]*8+[(-3,-2)],
                                     seed=seed,maxiter=220,popsize=12,polish=True,tol=1e-8)
        runs.append(dict(seed=seed,evaluations=int(result.nfev),saving=float(result.fun),
                         perturbation=result.x[:8].tolist(),threshold=float(10**result.x[8]),
                         parameters=perturbed(base,result.x[:8]).tolist(),
                         optimizer_success=bool(result.success)))
        print(seed,result.fun,'threshold',10**result.x[8],flush=True)
    record=dict(scope="continuous 8-dimensional +/-5% parameter box and KL threshold band",
                method="three seeded differential-evolution searches with local polish",
                runs=runs,elapsed_seconds=time.monotonic()-start,
                smallest_saving_found=min(r['saving'] for r in runs),
                establishes_global_minimum=False,certified_neighborhood=False)
    (ROOT/'physical_rate_adversarial_robustness.json').write_text(json.dumps(record,indent=2)+'\n')
