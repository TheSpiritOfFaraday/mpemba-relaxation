"""Search finite-time heating in TV, Hellinger and convex hinge divergences.
Every sampled generator is filtered by all eight original F/G inequalities.
Negative candidates are discovery only and require independent verification.
"""
import argparse,itertools,json,time,os
from pathlib import Path
import numpy as np
from cone_search import pairs,constraints
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=600)
args=parser.parse_args();P=Path(__file__).resolve().parent;start=time.time();last=start
rng=np.random.default_rng(111111);count=0;systems=0;best={};hits={}
for batch in itertools.count():
 if time.time()-start>args.seconds:break
 E=np.r_[0,np.cumsum(np.exp(rng.uniform(-5,1.7,3)))];q=np.exp(-E);q/=sum(q)
 K=np.zeros((500,4,4))
 for j in range(1,4):K[:,:j,j]=np.sort(np.exp(rng.uniform(-8,8,(500,j))),axis=1);K[:,j,:j]=K[:,:j,j]
 aa=np.array([q[i]*K[:,i,j] for i,j in pairs]).T
 good=(aa@constraints(q).T<=-1e-12*np.max(aa,axis=1)[:,None]).all(axis=1)
 K=K[good];nn=len(K)
 if not nn:continue
 L=K*q[None,None,:]
 for i in range(4):L[:,i,i]=-L[:,i,:].sum(axis=1)
 val,V=np.linalg.eigh(L*np.sqrt(q[None,:,None]/q[None,None,:]));VV=V[:,:,:3];ev=val[:,:3]
 systems+=nn
 beta=1+np.geomspace(.0001,100,20);p0=np.exp(-beta[:,None]*E);p0/=p0.sum(axis=1)[:,None]
 dp0=-p0*(E-(p0@E)[:,None]);co=np.einsum('sij,bi->sbj',VV,(p0-q)/np.sqrt(q));dc=np.einsum('sij,bi->sbj',VV,dp0/np.sqrt(q))
 for scaled in np.geomspace(.0001,5,25):
  tt=scaled/(-val[:,-2]);dec=np.exp(ev*tt[:,None]);delta=np.einsum('sij,sbj,sj->sbi',VV,co,dec)*np.sqrt(q)
  pt=q+delta;dt=np.einsum('sij,sbj,sj->sbi',VV,dc,dec)*np.sqrt(q)
  valid=(pt>1e-12).all(axis=2);scale=abs(dt).sum(axis=2);valid &=scale>1e-13
  tv=(dt*np.sign(delta)).sum(axis=2)/2
  hell=(dt*(1-np.sqrt(q/np.maximum(pt,1e-300)))).sum(axis=2)
  order=np.argsort(-pt/q,axis=2);ds=np.take_along_axis(dt,order,axis=2);cs=np.cumsum(ds,axis=2)[:,:,:3]
  hinge=cs.min(axis=2)
  count+=nn*len(beta)
  for name,deriv in [('TV',tv),('Hellinger',hell),('hinge',hinge)]:
   score=np.where(valid,deriv/np.maximum(scale,1e-300),np.inf);ix=np.unravel_index(np.argmin(score),score.shape);value=float(score[ix])
   if value<best.get(name,1):
    best[name]=value;i,b=ix
    data=dict(distance=name,E=E.tolist(),q=q.tolist(),K=K[i].tolist(),beta=float(beta[b]),time=float(tt[i]),scaled_time=float(scaled),score=value,derivative=float(deriv[ix]),p=pt[i,b].tolist(),dp=dt[i,b].tolist(),ratio_order=order[i,b].tolist(),cumulative_derivatives=cs[i,b].tolist(),eigenvalues=val[i].tolist(),seed=111111,evaluations=count,status='numerical discovery')
    (P/f'transient_{name}_candidate.json').write_text(json.dumps(data,indent=2)+'\n')
    if value<-1e-6:
     hits[name]=data;print('CANDIDATE',name,json.dumps(data),flush=True)
  if len(hits)==3:break
 if len(hits)==3:break
 if time.time()-last>30:
  status=dict(pid=os.getpid(),status='running',systems=systems,evaluations=count,best=best,hits=list(hits),elapsed_seconds=time.time()-start)
  (P/'transient_distances_status.json').write_text(json.dumps(status,indent=2)+'\n');print('PROGRESS',json.dumps(status),flush=True);last=time.time()
status=dict(pid=os.getpid(),status='finished',systems=systems,evaluations=count,best=best,hits=list(hits),elapsed_seconds=time.time()-start)
(P/'transient_distances_status.json').write_text(json.dumps(status,indent=2)+'\n');print('FINISHED',json.dumps(status),flush=True)
