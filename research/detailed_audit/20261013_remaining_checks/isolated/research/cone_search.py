"""Reproducible ray-mixture search, including signed directions plus reset.
Seeds, ranges and counts recorded by output. This is discovery, not certification.
"""
import numpy as np, itertools,json
from pathlib import Path
pairs=list(itertools.combinations(range(4),2)); trips=list(itertools.combinations(range(4),3))
def constraints(q):
 rows=[]
 for i,j,k in trips:
  t=q[j]/q[i];u=q[k]/q[j]
  f=np.zeros(6);g=f.copy()
  f[pairs.index((i,k))]=1+t*u/2;f[pairs.index((j,k))]=-1-u/2;f[pairs.index((i,j))]=-(1-t)
  g[pairs.index((i,k))]=t;g[pairs.index((j,k))]=-1
  rows.extend([f,g])
 return np.array(rows)
def rays(q):
 C=constraints(q)[:,1:];out=[]
 for ids in itertools.combinations(range(8),4):
  U,s,V=np.linalg.svd(C[list(ids)])
  if s[-1]<1e-9:continue
  for v in [V[-1],-V[-1]]:
   if (C@v).max()<1e-8 and not any(np.linalg.norm(v-r)<1e-6 for r in out):out.append(v)
 return np.pad(np.array(out),((0,0),(1,0)))
def generator(q,a):
 M=np.zeros((4,4))
 for (i,j),r in zip(pairs,a): M[i,j]=r;M[j,i]=r*q[j]/q[i]
 M-=np.diag(M.sum(axis=0));return M

def run(seed=10411,energy_cases=300,mixtures=1000):
 rng=np.random.default_rng(seed);count=0;best=-1;candidate=None
 betas=np.linspace(.001,.999,150)
 for z in range(energy_cases):
  gaps=np.exp(rng.uniform(-4,2,3));E=np.r_[0,np.cumsum(gaps)];q=np.exp(-E);rs=rays(q)
  if len(rs)==0:continue
  weights=np.exp(rng.uniform(-9,3,(mixtures,len(rs))))
  aa=weights@rs
  # Signed rates are safe for eigenvector discovery; add reset to restore positivity.
  reset=np.array([q[i] for i,j in pairs]);rho=np.maximum(0,(-aa/reset).max(axis=1))+1
  aa+=rho[:,None]*reset
  M=np.array([generator(q,a) for a in aa]);S=M*np.sqrt(q[None,None,:]/q[None,:,None])
  lam,V=np.linalg.eigh(S);w=V[:,:,-2]/np.sqrt(q)
  probs=np.exp(-betas[:,None]*E);probs/=probs.sum(axis=1)[:,None]
  b=w@probs.T;der=w@(-probs*(E-(probs@E)[:,None])).T
  score=b*der # positive => cooling anomaly locally
  idx=np.unravel_index(np.argmax(score),score.shape);val=score[idx];count+=mixtures
  if val>best:
   best=float(val);a=aa[idx[0]]
   candidate=dict(E=E.tolist(),q=q.tolist(),a=a.tolist(),w=w[idx[0]].tolist(),eigenvalues=lam[idx[0]].tolist(),beta=float(betas[idx[1]]),score=best,seed=seed,cases=count)
  if val>1e-7:
   Path('research/cone_candidate.json').write_text(json.dumps(candidate,indent=2));print(json.dumps(candidate),flush=True);return
  if z%20==0:print('progress',z,count,best,flush=True)
 Path('research/cone_search_result.json').write_text(json.dumps(dict(seed=seed,count=count,best=best,candidate=candidate),indent=2))
if __name__=='__main__':run()
