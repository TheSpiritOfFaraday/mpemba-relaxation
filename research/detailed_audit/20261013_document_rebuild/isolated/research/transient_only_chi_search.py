"""Find two finite-time chi-squared crossings with normal eventual ordering.
Full original F/G constraints enforced numerically; discoveries need certificates.
"""
import numpy as np,itertools,time,json,os,argparse
from pathlib import Path
from cone_search import pairs,constraints
p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=600);args=p.parse_args()
rng=np.random.default_rng(111120);start=time.time();count=0;best=0.;P=Path(__file__).resolve().parent
for batch in itertools.count():
 if time.time()-start>args.seconds:break
 E=np.r_[0,np.cumsum(np.exp(rng.uniform(-6,2,3)))];q=np.exp(-E);q/=sum(q)
 K=np.zeros((1000,4,4))
 for j in range(1,4):K[:,:j,j]=np.sort(np.exp(rng.uniform(-8,8,(1000,j))),axis=1);K[:,j,:j]=K[:,:j,j]
 aa=np.array([q[i]*K[:,i,j] for i,j in pairs]).T
 good=(aa@constraints(q).T<=-1e-11*np.max(aa,axis=1)[:,None]).all(axis=1)
 aa=aa[good];K=K[good];nn=len(K)
 if not nn:continue
 L=K*q[None,None,:]
 for i in range(4):L[:,i,i]=-L[:,i,:].sum(axis=1)
 val,V=np.linalg.eigh(L*np.sqrt(q[None,:,None]/q[None,None,:]));val=val[:,:3][:,::-1];V=V[:,:,:3][:,:,::-1]
 for bh,bc in [(.001,.1),(.001,.5),(.1,.5),(1/3,2/3),(.001,.8),(.01,.05)]:
  pp=np.exp(-np.array([bh,bc])[:,None]*E);pp/=pp.sum(axis=1)[:,None]
  coef=np.einsum('sij,bi->sbj',V,(pp-q)/np.sqrt(q));abc=coef[:,0,:]**2-coef[:,1,:]**2
  A,B,C=abc.T;a=2*(val[:,0]-val[:,1]);b=2*(val[:,0]-val[:,2])
  valid=(A>1e-8*np.sum(abs(abc),axis=1)) & (B<0) & (C>0) & (b>a+1e-8*np.maximum(b,1e-300))
  ratio=np.where(valid,-b*C/np.minimum(a*B,-1e-300),1);valid &=ratio>1
  tt=np.log(np.maximum(ratio,1))/(b-a+1e-300)
  F=A+B*np.exp(-a*tt)+C*np.exp(-b*tt);score=F/np.maximum(np.sum(abs(abc),axis=1),1e-300)
  score[~valid]=1;ix=np.argmin(score);count+=nn
  if score[ix]<best:
   best=float(score[ix]);data=dict(E=E.tolist(),q=q.tolist(),rates=aa[ix].tolist(),beta_hot=bh,beta_cold=bc,coefficients=abc[ix].tolist(),eigenvalues=val[ix].tolist(),critical_time=float(tt[ix]),score=best,seed=111120,evaluations=count,status='numerical candidate')
   (P/'transient_only_chi_candidate.json').write_text(json.dumps(data,indent=2)+'\n');print('BEST',json.dumps(data),flush=True)
  if best<-1e-6:break
 if best<-1e-6:break
 if batch%100==0:print('PROGRESS',batch,count,best,flush=True)
out=dict(pid=os.getpid(),status='finished',evaluations=count,best=best,elapsed_seconds=time.time()-start)
(P/'transient_only_chi_status.json').write_text(json.dumps(out,indent=2)+'\n');print('FINISHED',json.dumps(out),flush=True)
