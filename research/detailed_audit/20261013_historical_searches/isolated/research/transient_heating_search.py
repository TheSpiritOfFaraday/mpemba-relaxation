"""Adversarial finite-time KL heating search under G-only ordering.
This tests a possible extension NOT claimed by the contrast theorem.
Numerical candidates require separate high-precision or exact verification.
"""
import argparse,itertools,json,time,os
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=1200);args=parser.parse_args()
rng=np.random.default_rng(111032);start=time.time();last=start;count=0;best=1.;bestdata=None;hits=0
P=Path(__file__).resolve().parent
for batch in itertools.count():
 if time.time()-start>args.seconds:break
 n=[4,5,6,8][batch%4];E=np.r_[0,np.cumsum(10**rng.uniform(-2.5,.5,n-1))]
 q=np.exp(-E);q/=sum(q)
 if min(q)<1e-8:continue
 nn=200;K=np.zeros((nn,n,n))
 for j in range(1,n):
    K[:,:j,j]=np.sort(np.exp(rng.uniform(-7,7,(nn,j))),axis=1);K[:,j,:j]=K[:,:j,j]
 L=K*q[None,None,:]
 for i in range(n):L[:,i,i]=-L[:,i,:].sum(axis=1)
 S=L*np.sqrt(q[None,:,None]/q[None,None,:]);vals,V=np.linalg.eigh(S)
 if np.any(vals[:,-2]>=-1e-10):continue
 # Only nonstationary modes; exact bath supplies the stationary term.
 VV=V[:,:,:-1];ev=vals[:,:-1]
 bs=1+np.r_[np.geomspace(.0001,1,8),np.geomspace(1.5,50,8)]
 pp=np.exp(-bs[:,None]*E);pp/=pp.sum(axis=1)[:,None]
 dpp=-pp*(E-(pp@E)[:,None])
 co=np.einsum('sij,b i->sbj',VV,(pp-q)/np.sqrt(q))
 dc=np.einsum('sij,b i->sbj',VV,dpp/np.sqrt(q))
 for scaled in np.geomspace(.0001,6,24):
    tt=scaled/(-vals[:,-2]);dec=np.exp(ev*tt[:,None])
    delta=np.einsum('sij,sbj,sj->sbi',VV,co,dec)*np.sqrt(q)
    pt=q+delta;dt=np.einsum('sij,sbj,sj->sbi',VV,dc,dec)*np.sqrt(q)
    valid=(pt>0).all(axis=2)
    logs=np.log(np.maximum(pt,1e-300)/q);terms=dt*logs
    derivative=terms.sum(axis=2);scale=abs(terms).sum(axis=2)
    score=np.where(valid & (scale>1e-14),derivative/np.maximum(scale,1e-300),np.inf)
    ix=np.unravel_index(np.argmin(score),score.shape);value=float(score[ix]);count+=nn*len(bs)
    if value<best:
        best=value;i,b=ix;bestdata=dict(n=n,E=E.tolist(),q=q.tolist(),K=K[i].tolist(),beta=float(bs[b]),time=float(tt[i]),scaled_time=float(scaled),score=best,derivative=float(derivative[ix]),p=pt[i,b].tolist(),eigenvalues=vals[i].tolist())
    if value<-1e-7:
        hits+=1;(P/f'transient_heating_candidate_{hits}.json').write_text(json.dumps(bestdata,indent=2)+'\n')
        print('CANDIDATE',hits,json.dumps(bestdata),flush=True)
        if hits>=5:break
 if hits>=5:break
 if time.time()-last>30:
    status=dict(pid=os.getpid(),seed=111032,elapsed_seconds=time.time()-start,evaluations=count,best_score=best,hits=hits,best=bestdata,status='running')
    (P/'transient_heating_status.json').write_text(json.dumps(status,indent=2)+'\n');print('PROGRESS',count,best,hits,flush=True);last=time.time()
status=dict(pid=os.getpid(),seed=111032,elapsed_seconds=time.time()-start,evaluations=count,best_score=best,hits=hits,best=bestdata,status='finished')
(P/'transient_heating_status.json').write_text(json.dumps(status,indent=2)+'\n');print('FINISHED',json.dumps(status),flush=True)
