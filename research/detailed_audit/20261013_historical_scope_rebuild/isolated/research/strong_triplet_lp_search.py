"""Adversarial strong-effect search with prescribed thermal zero and PSD cuts.
Numerical discovery only. Reversible four-state generators, all original F/G<=0,
M^T w=-w, and the other two nonstationary rates required to exceed 1.01.
The cutting-plane LP imposes 2x2 semidefinite positivity on the fast subspace.
All candidate claims need separate exact certification.
"""
import argparse,itertools,json,time,os
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from cone_search import constraints,generator
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=10800);parser.add_argument('--seed',type=int,default=111026);args=parser.parse_args()
rng=np.random.default_rng(args.seed);start=time.time();count=lpcount=feasible=0;bestgap=-1e9;best=None;last=start
for count in itertools.count(1):
 if time.time()-start>args.seconds:break
 if count%4==0:E=np.array([0,np.log(200/199),np.log(2.5),np.log(400)])
 else:E=np.r_[0,np.cumsum(10**rng.uniform(-3,.9,3))]
 q=np.exp(-E);pi=q/sum(q)
 if pi.min()<1e-10:continue
 if count%2:beta=10**rng.uniform(-3,-.005)
 else:beta=1+10**rng.uniform(-2,1.5)
 pb=np.exp(-beta*E);pb/=sum(pb)
 _,ss,V=np.linalg.svd(np.array([pi,pb]),full_matrices=True)
 if ss[-1]<1e-9:continue
 theta=rng.uniform(0,2*np.pi);w=V[2]*np.cos(theta)+V[3]*np.sin(theta);w/=max(abs(w))
 sq=np.sqrt(pi);v=sq*w;v/=np.linalg.norm(v)
 _,_,VV=np.linalg.svd(np.array([sq,v]),full_matrices=True);U=VV[2:].T
 Ms=np.array([generator(q,np.eye(6)[i]) for i in range(6)])
 eq=np.array([m.T@w for m in Ms]).T
 Cs=constraints(q)
 Ss=Ms*sq[None,None,:]/sq[None,:,None]
 fast=np.array([-U.T@S@U for S in Ss])
 # Directions of the semidefinite cone; add violating eigenvectors adaptively.
 dirs=[np.array([1.,0]),np.array([0,1.]),np.array([1.,1.])/np.sqrt(2),np.array([1.,-1.])/np.sqrt(2)]
 result=None
 for cut in range(12):
  extra=np.array([-np.einsum('i,kij,j->k',d,fast,d) for d in dirs])
  A=np.r_[Cs,extra];b=np.r_[np.full(8,-1e-9),np.full(len(dirs),-1.01)]
  res=linprog(np.ones(6),A_ub=A,b_ub=b,A_eq=eq,b_eq=-w,bounds=[(1e-8,1e7)]*6,method='highs')
  lpcount+=1
  if not res.success:break
  fastmat=np.einsum('k,kij->ij',res.x,fast);ev,vec=np.linalg.eigh(fastmat);gap=float(ev[0]-1)
  if gap>bestgap:
   bestgap=gap;best=dict(E=E.tolist(),beta=beta,w=w.tolist(),rates=res.x.tolist(),fast_eigenvalues=ev.tolist(),triplet_residuals=(Cs@res.x).tolist(),eq_residual=float(max(abs(eq@res.x+w))),count=count)
  if ev[0]>=1.01-1e-8:
   result=res;break
  dirs.append(vec[:,0])
 if result is not None:
  rates=result.x;M=generator(q,rates);ev=np.linalg.eigvalsh(M*sq[None,:]/sq[:,None]);err=max(abs(M.T@w+w));fg=Cs@rates
  if ev[2]<-.99999 and ev[1]<-1.005 and err<1e-7 and max(fg)<1e-7:
   feasible+=1;out=dict(best or {},seed=args.seed,elapsed=time.time()-start,eigenvalues=ev.tolist(),thermal_zero_residual=float(w@pb),bath_zero_residual=float(w@pi),label='NUMERICAL CANDIDATE ONLY')
   (P/f'strong_triplet_lp_candidate_{feasible}.json').write_text(json.dumps(out,indent=2)+'\n')
   print('CANDIDATE',feasible,json.dumps(out),flush=True)
   if feasible>=5:break
 if time.time()-last>30:
  status=dict(pid=os.getpid(),seed=args.seed,elapsed_seconds=time.time()-start,systems=count,lp_calls=lpcount,candidates=feasible,best_fast_gap=bestgap,best=best,status='running')
  (P/'strong_triplet_lp_status.json').write_text(json.dumps(status,indent=2)+'\n');print('PROGRESS',count,lpcount,feasible,bestgap,flush=True);last=time.time()
status=dict(pid=os.getpid(),seed=args.seed,elapsed_seconds=time.time()-start,systems=count,lp_calls=lpcount,candidates=feasible,best_fast_gap=bestgap,best=best,status='finished')
(P/'strong_triplet_lp_status.json').write_text(json.dumps(status,indent=2)+'\n');print('FINISHED',json.dumps(status),flush=True)
