"""Search for an anomalous cooling interval with a positive left endpoint.
The linear program first tests the larger positive contrast-mode cone.
Any hit is then tested for full-triplet inverse realizability. Discovery only.
"""
import numpy as np
from scipy.optimize import linprog
from mode_realizability import contrast_basis,realize
from pathlib import Path
import json,time
rng=np.random.default_rng(111107);start=time.time();hits=[];tests=0
for it in range(2000):
 E=np.r_[0,np.cumsum(np.exp(rng.uniform(-6,4.5,3)))];q=np.exp(-E)
 if E[-1]>250:continue
 H=contrast_basis(q)
 # Normalize each basis function by its own initial amplitude. Without this,
 # large energy spans make the LP ignore the smaller contrast coordinates.
 H/=-(H.T@np.ones(4)/4)[None,:]
 def deriv(b):
  p=np.exp(-b*E);p/=sum(p);return H.T@(-p*(E-p@E))
 d0=deriv(0);d0/=np.linalg.norm(d0)
 for beta in [.005,.02,.07,.2,.5,.8]:
  db=deriv(beta);db/=np.linalg.norm(db)
  res=linprog(db,A_ub=[-d0],b_ub=[-.001],A_eq=[[1,1,1]],b_eq=[1],bounds=[(0,None)]*3,method='highs');tests+=1
  if res.success and res.fun<-.0001:
   w=H@res.x;r=realize(q,w)
   hits.append(dict(E=E.tolist(),w=w.tolist(),beta=beta,objective=float(res.fun),coefficients=res.x.tolist(),realizable=r['feasible']))
   if len(hits)>=10:break
 if len(hits)>=10:break
out=dict(seed=111107,systems=it+1,lp_tests=tests,seconds=time.time()-start,hits=hits,scope='Exploration of cooling intervals with positive left endpoint; not proof')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
