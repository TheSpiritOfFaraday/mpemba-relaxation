"""Adversarial LP cone search for general-n first-three arithmetic inequality."""
import argparse,json,time,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

def setup(q):
 n=len(q);pairs=list(itertools.combinations(range(n),2));idx={p:i for i,p in enumerate(pairs)};rows=[]
 for i,j,k in itertools.combinations(range(n),3):
  f=np.zeros(len(pairs));g=f.copy()
  f[idx[i,k]]=q[i]+q[k]/2;f[idx[j,k]]=-q[j]-q[k]/2;f[idx[i,j]]=-(q[i]-q[j])
  g[idx[i,k]]=1;g[idx[j,k]]=-1
  rows.extend([f/np.max(np.abs(f)),g])
 rows=np.array(rows)
 Z=np.zeros((n,n-1))
 for j in range(n-1):
  Q=np.sum(q[:j+1]);Qn=Q+q[j+1]
  Z[:j+1,j]=np.sqrt(q[:j+1]*q[j+1]/(Q*Qn));Z[j+1,j]=-np.sqrt(Q/Qn)
 edges=np.array([np.outer(np.sqrt(q[j])*Z[i]-np.sqrt(q[i])*Z[j],np.sqrt(q[j])*Z[i]-np.sqrt(q[i])*Z[j]) for i,j in pairs])
 d=(q[0]-q[1])*np.sqrt(q[2])/(2*np.sqrt(q[0]*q[1]*sum(q[:3])))
 return pairs,rows,Z,edges,d

def evaluate(k,q,data):
 pairs,rows,Z,edges,d=data
 A=-np.tensordot(k,edges,axes=(0,0));eig,u=np.linalg.eigh(A);v=u[:,-1]
 if sum(v)<0:v=-v
 ratio=(v[1]-d*v[0])/(abs(v[1])+d*abs(v[0])+1e-300)
 return ratio,eig[-1]-eig[-2],v,eig,A

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=int,default=180);ap.add_argument('--seed',type=int,default=917352);args=ap.parse_args();rng=np.random.default_rng(args.seed);start=time.time();counts=dict(q=0,lp=0,lp_success=0,modes=0,feasible_modes=0);best=None
 while time.time()-start<args.seconds:
  n=5+counts['q']%2;regime=counts['q']%6
  gaps=10**rng.uniform(-3,1,n-1)
  if regime==0:gaps[:2]=10**rng.uniform(-4,-1,2)
  if regime==1:gaps[2:]=10**rng.uniform(-4,-1,n-3)
  if regime==2:gaps[2]=rng.uniform(10,40)
  if regime==3:gaps[0]=rng.uniform(5,20)
  if regime==4:gaps[:]=rng.uniform(.01,3,n-1)
  if regime==5:gaps[:]=10**rng.uniform(-3,.8,n-1)
  E=np.r_[0,np.cumsum(gaps)];q=np.exp(-E);q/=sum(q);data=setup(q);pairs,rows,Z,edges,d=data;vertices=[];counts['q']+=1
  bounds=[(0,0)]+[(-1,1)]*(len(pairs)-1)
  for _ in range(20):
   c=rng.normal(size=len(pairs));sol=linprog(c,A_ub=rows,b_ub=np.zeros(len(rows)),bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9});counts['lp']+=1
   if not sol.success:continue
   counts['lp_success']+=1
   if np.max(rows@sol.x)>2e-9:continue
   vertices.append(sol.x)
  if not vertices:continue
  vertices=np.array(vertices)
  for t in range(100):
   if t<len(vertices):k=vertices[t]
   else:
    # Broad mixtures explore both interior and nearly extreme directions.
    coeff=10**rng.uniform(-8,0,len(vertices));k=coeff@vertices/sum(coeff)
   ratio,gap,v,eig,A=evaluate(k,q,data);counts['modes']+=1
   if gap<1e-9 or min(v)<-1e-6 or abs(v[1])+d*abs(v[0])<1e-7:continue
   counts['feasible_modes']+=1
   if best is None or ratio<best['ratio']:
    best=dict(n=n,E=E.tolist(),q=q.tolist(),reduced_rates=k.tolist(),ratio=float(ratio),gap=float(gap),perron=v.tolist(),eigenvalues=eig.tolist(),max_constraint_violation=float(np.max(rows@k)),regime=regime)
   if ratio < -1e-5 and gap>1e-7 and v[1]-d*v[0] < -1e-8:
    out=dict(candidate=True,seed=args.seed,elapsed=time.time()-start,counts=counts,best=best)
    Path('research/general_n_arithmetic_candidate.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True);return
  if counts['q']%100==0:print(json.dumps(dict(elapsed=time.time()-start,counts=counts,best_ratio=best['ratio'] if best else None)),flush=True)
 out=dict(candidate=False,seed=args.seed,seconds=args.seconds,elapsed=time.time()-start,counts=counts,best=best)
 Path('research/general_n_arithmetic_search_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
if __name__=='__main__':main()
