"""Adversarial-corner refinement within the unchanged declared physical model."""
import itertools
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
import search_physical_rate_model as search
from physical_rate_model import Evolution

ROOT=Path(__file__).resolve().parent
corners=np.array(list(itertools.product([-1,1],repeat=8)))


def scores(x):
    return np.array([search.score(search.perturb(x,1+.05*s)) for s in corners])


if __name__=='__main__':
    start=time.monotonic()
    x=np.array(json.loads((ROOT/'physical_rate_search_robust_screen.json').read_text())['result']['parameters'])
    bounds=[(.11,2),(.11,2),(.11,4.5),(2.11,5.7),(2.11,5.7),(2.11,5.7),(1.4,5),(2,9.5)]
    stages=[]
    active=set()
    for stage in range(3):
        before=scores(x)
        active.update(np.argsort(before)[-12:].tolist())
        def objective(y):
            if sum(y[:3])>4.76 or y[7]<y[6]+.3:
                return 10+max(sum(y[:3])-4.76,y[6]+.3-y[7])
            r=search.box_rate_bound(y)
            if r>995:return 10+np.log(r/995)
            return max(search.score(search.perturb(y,1+.05*corners[i])) for i in sorted(active))
        result=differential_evolution(objective,bounds,seed=12102028+stage,
                                     maxiter=120,popsize=10,polish=True,x0=x,tol=1e-6)
        proposed=result.x
        after=scores(proposed)
        # Keep the best candidate on all corners, not just the active training set.
        accepted=bool(max(after)<max(before))
        if accepted:x=proposed
        stages.append(dict(stage=stage,seed=12102028+stage,
                           active_corner_indices=sorted(active),
                           evaluations=int(result.nfev),score=float(result.fun),
                           all_corner_score=float(max(after)),accepted=accepted,
                           parameters=proposed.tolist()))
        active.update(np.argsort(after)[-12:].tolist())
        receipt=dict(model="two-basin shared barrier excess",stages=stages,
                     elapsed_seconds=time.monotonic()-start,
                     all_corner_score=float(max(scores(x))),
                     result=Evolution(x).metrics(),certified_neighborhood=False)
        (ROOT/'physical_rate_refined.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print('stage',stage,'worst corner residual',receipt['all_corner_score'],
              'parameters',x.tolist(),flush=True)
