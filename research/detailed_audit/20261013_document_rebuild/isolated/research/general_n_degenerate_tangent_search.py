"""Target splitting of repeated slow eigenspace at one-fast adjacent-pair boundary."""
import json,time,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
from general_n_arithmetic_search import setup,evaluate
rng=np.random.default_rng(914728);start=time.time();counts=dict(baths=0,boundaries=0,lp=0,modes=0,retained=0);best=1
while time.time()-start<120:
 n=int(rng.choice([5,6]));gaps=10**rng.uniform(-3,.7,n-1);E=np.r_[0,np.cumsum(gaps)];q=np.exp(-E);q/=sum(q);data=setup(q);pairs,rows,Z,edges,d=data;counts['baths']+=1
 for pair in range(n-1):
  core=np.zeros(len(pairs));core[pairs.index((pair,pair+1))]=1/(q[pair]+q[pair+1]);active=rows[np.abs(rows@core)<1e-10];f=np.sqrt(q[pair+1]/(q[pair]+q[pair+1]))*Z[pair]-np.sqrt(q[pair]/(q[pair]+q[pair+1]))*Z[pair+1];V=null_space(f.reshape(1,-1));counts['boundaries']+=1;vertices=[]
  for _ in range(10):
   sol=linprog(rng.normal(size=len(pairs)),A_ub=active,b_ub=np.zeros(len(active)),bounds=[(-1,1)]*len(pairs),method='highs',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10});counts['lp']+=1
   if sol.success and np.max(active@sol.x)<1e-9:vertices.append(sol.x)
  if not vertices:continue
  for _ in range(40):
   weights=10**rng.uniform(-5,0,len(vertices));k=weights@np.array(vertices)/sum(weights);A=-np.tensordot(k,edges,axes=(0,0));eff=V.T@A@V;ev,v=np.linalg.eigh(eff);u=V@v[:,-1]
   if sum(u)<0:u=-u
   counts['modes']+=1
   if ev[-1]-ev[-2]<1e-7 or min(u)<-1e-8 or abs(u[1])+d*abs(u[0])<1e-7:continue
   counts['retained']+=1;ratio=(u[1]-d*u[0])/(abs(u[1])+d*abs(u[0]));best=min(best,ratio)
   if ratio< -1e-5 and u[1]-d*u[0]<-1e-8:
    out=dict(candidate=True,counts=counts,n=n,q=q.tolist(),E=E.tolist(),fast_pair=[pair+1,pair+2],perturb_reduced_rates=k.tolist(),perron=u.tolist(),ratio=float(ratio),effective_gap=float(ev[-1]-ev[-2]),elapsed=time.time()-start)
    Path('research/general_n_degenerate_tangent_candidate.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));raise SystemExit
out=dict(candidate=False,counts=counts,best_relative_margin=float(best),seed=914728,elapsed=time.time()-start,scope='Leading simple-slow branch obtained by splitting repeated slow boundary of reset plus one adjacent fast pair; perturbations obey all active original F/G inequalities.')
Path('research/general_n_degenerate_tangent_search_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
