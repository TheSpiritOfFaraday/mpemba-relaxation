"""Numerical adversarial fixed-mode LP; no infeasibility claim constitutes proof."""
import numpy as np
from scipy.optimize import linprog
from itertools import combinations
import json,time
rng=np.random.default_rng(723401)
start=time.time()
def test(q,x):
 n=len(q);m=n-1;H=np.zeros((n,m));Q=np.cumsum(q)
 for j in range(m):H[:j+1,j]=1;H[j+1,j]=-Q[j]/q[j+1]
 eta=np.diag(H.T@np.diag(q)@H);edges=list(combinations(range(m),2)); Ms=[]
 for i,j in edges:
  K=np.zeros((m,m));K[i,j]=1/eta[i];K[j,i]=1/eta[j]
  K[i,i]=-x[j]/(eta[i]*x[i]);K[j,j]=-x[i]/(eta[j]*x[j])
  Ms.append((H@K@np.diag(1/eta)@H.T@np.diag(q)).T)
 Ms=np.array(Ms);rows=[]
 for i,j,k in combinations(range(n),3):
  rows.append(Ms[:,i,k]*(1+q[k]/(2*q[i]))-Ms[:,j,k]*(1+q[k]/(2*q[j]))-Ms[:,i,j]*(1-q[j]/q[i]))
  rows.append(Ms[:,i,k]*q[j]/q[i]-Ms[:,j,k])
 A=np.array(rows);scale=np.max(abs(A),axis=1);A=A/np.maximum(scale,1e-100)[:,None]
 # maximize a shared positive margin on all F/G and on every graph edge
 N=len(edges);ub=np.r_[np.c_[A,np.ones(len(A))],np.c_[-np.eye(N),np.ones(N)]]
 sol=linprog(np.r_[np.zeros(N),-1.],A_ub=ub,b_ub=np.zeros(len(ub)),A_eq=[np.r_[np.ones(N),0]],b_eq=[1],bounds=[(0,None)]*N+[(None,None)],method='highs',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
 if sol.success and sol.x[-1]>1e-8 and np.min(sol.x[:-1])>sol.x[-1]*.9 and np.max(A@sol.x[:-1]) < -sol.x[-1]*.9:
  return dict(q=q.tolist(),x=x.tolist(),edges=edges,ell=sol.x[:-1].tolist(),margin=sol.x[-1],M_shift=np.einsum('i,ijk->jk',sol.x[:-1],Ms).tolist())
 return None
for z in range(30000):
 n=5 if z<20000 else 6
 if z%3==0: gaps=rng.uniform(.02,1.5,n-1)
 elif z%3==1:gaps=10**rng.uniform(-3,.4,n-1)
 else:gaps=rng.uniform(.1,1,n-1)
 q=np.r_[1,np.exp(-np.cumsum(gaps))]
 s=(q[0]-q[1])*q[2]/(2*q[1]*sum(q[:3]))
 x=10**rng.uniform(-4,1,n-1);x[0]=1;x[1]=s*rng.uniform(.05,.999)
 # concentrate tail modes also below the tiny middle coordinate
 if z%2:x[2:]*=s
 out=test(q,x)
 if out:
  out.update(iteration=z,seed=723401,elapsed=time.time()-start)
  open('research/higher_dimension_first3_candidate.json','w').write(json.dumps(out,indent=2));print('FOUND',z,out['margin'],flush=True);break
 if z%1000==0:print(z,round(time.time()-start,1),flush=True)
else:print('NO_STRICT_COUNTEREXAMPLE',z+1,round(time.time()-start,1),flush=True)
