"""Test F123 after removing minimal reset and tracing all higher states."""
import itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from general_n_arithmetic_search import setup
rng=np.random.default_rng(919536);start=time.time();count=0;best=-1e99
for trial in range(150000):
 n=5+trial%2;E=np.r_[0,np.cumsum(10**rng.uniform(-3,1,n-1))];q=np.exp(-E);q/=sum(q);pairs,rows,*_=setup(q)
 sols=[]
 for _ in range(3):
  sol=linprog(rng.normal(size=len(pairs)),A_ub=rows,b_ub=np.zeros(len(rows)),bounds=[(0,0)]+[(-1,1)]*(len(pairs)-1),method='highs')
  if sol.success and max(rows@sol.x)<1e-9:sols.append(sol.x)
 if not sols:continue
 coeff=10**rng.uniform(-5,0,len(sols));k=coeff@np.array(sols)/sum(coeff);k-=min(k);K=np.zeros((n,n))
 for (i,j),z in zip(pairs,k):K[i,j]=K[j,i]=z
 M=q[:,None]*K;np.fill_diagonal(M,-np.sum(M,axis=0))
 if np.linalg.cond(M[3:,3:])>1e12:continue
 T=M[:3,:3]-M[:3,3:]@np.linalg.solve(M[3:,3:],M[3:,:3]);Kt=T/q[:3,None]
 F=(q[0]+q[2]/2)*Kt[0,2]-(q[1]+q[2]/2)*Kt[1,2]-(q[0]-q[1])*Kt[0,1];count+=1
 if F>best:best=F
 if F>1e-8:
  out=dict(candidate=True,seed=919536,trial=trial,count=count,elapsed=time.time()-start,n=n,q=q.tolist(),K=K.tolist(),trace_K=Kt.tolist(),trace_F=float(F),original_max_constraint=float(max(rows@k)))
  Path('research/general_n_reset_trace_candidate.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));raise SystemExit
out=dict(candidate=False,count=count,best=float(best),elapsed=time.time()-start)
Path('research/general_n_reset_trace_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
