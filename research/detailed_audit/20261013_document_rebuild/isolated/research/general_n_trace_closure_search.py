"""Search for original F closure failure under elimination of highest state."""
import itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from general_n_arithmetic_search import setup
rng=np.random.default_rng(871363);start=time.time();count=0;best=-1e99
for trial in range(30000):
 n=4+trial%3;E=np.r_[0,np.cumsum(10**rng.uniform(-3,1,n-1))];q=np.exp(-E);q/=sum(q);pairs,rows,*_=setup(q)
 sol=linprog(rng.normal(size=len(pairs)),A_ub=rows,b_ub=np.zeros(len(rows)),bounds=[(0,0)]+[(-1,1)]*(len(pairs)-1),method='highs')
 if not sol.success or max(rows@sol.x)>1e-8:continue
 k=sol.x-min(sol.x)+(0 if trial%2 else 10**rng.uniform(-3,2));K=np.zeros((n,n))
 for (i,j),z in zip(pairs,k):K[i,j]=K[j,i]=z
 d=q[:-1]@K[:-1,-1]
 if d<1e-12:continue
 Kt=K[:-1,:-1]+q[-1]*np.outer(K[:-1,-1],K[:-1,-1])/d
 for i,j,l in itertools.combinations(range(n-1),3):
  Ft=(q[i]+q[l]/2)*Kt[i,l]-(q[j]+q[l]/2)*Kt[j,l]-(q[i]-q[j])*Kt[i,j];count+=1
  if Ft>best:best=Ft
  if Ft>1e-7:
   out=dict(candidate=True,seed=871363,trial=trial,count=count,elapsed=time.time()-start,n=n,q=q.tolist(),K=K.tolist(),trace_K=Kt.tolist(),bad_triplet=[i+1,j+1,l+1],trace_F=float(Ft),min_k=float(min(k)),original_max_constraint=float(max(rows@k)))
   Path('research/general_n_trace_closure_candidate.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));raise SystemExit
out=dict(candidate=False,count=count,best=float(best),elapsed=time.time()-start)
Path('research/general_n_trace_closure_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
