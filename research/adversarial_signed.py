import numpy as np,json,time
from itertools import combinations
rng=np.random.default_rng(742908); pairs=list(combinations(range(4),2)); trips=list(combinations(range(4),3));start=time.time();accepted=count=nonmon=0
for batch in range(3000):
 n=10000
 E=np.column_stack([np.zeros(n),np.cumsum(10**rng.uniform(-3,.7,(n,3)),axis=1)])
 q=np.exp(-E);p=q/q.sum(axis=1,keepdims=True)
 rates=rng.normal(size=(n,6))*10**rng.uniform(-3,3,(n,6))
 if batch%2: rates*=q[:,[i for i,j in pairs]]
 keep=np.ones(n,bool)
 for i,j,k in trips:
  A=rates[:,pairs.index((i,k))];B=rates[:,pairs.index((j,k))];C=rates[:,pairs.index((i,j))];t=q[:,j]/q[:,i];u=q[:,k]/q[:,j]
  F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
  keep&=(F<=0)&(G<=0)
 E=E[keep];q=q[keep];p=p[keep];rates=rates[keep];m=len(E);count+=n;accepted+=m
 if not m:continue
 S=np.zeros((m,4,4));diag=np.zeros((m,4))
 for z,(i,j) in enumerate(pairs):
  a=rates[:,z];up=a*q[:,j]/q[:,i];S[:,i,j]=S[:,j,i]=a*np.sqrt(q[:,j]/q[:,i]);diag[:,j]-=a;diag[:,i]-=up
 for i in range(4):S[:,i,i]=diag[:,i]
 s=np.sqrt(p);v=s.copy();v[:,0]-=1;v/=np.linalg.norm(v,axis=1,keepdims=True)
 H=np.eye(4)[None,:,:]-2*v[:,:,None]*v[:,None,:];B=H[:,:,1:]
 T=np.swapaxes(B,1,2)@S@B
 vals,vec=np.linalg.eigh(T);w=np.einsum('nij,nj->ni',B,vec[:,:,2])/s
 dw=np.diff(w,axis=1);nmon=(dw.min(axis=1)<-1e-8)&(dw.max(axis=1)>1e-8);nonmon+=nmon.sum()
 if np.any(nmon) and batch==0: print('nonmon',json.dumps(dict(E=E[nmon][0].tolist(),rates=rates[nmon][0].tolist(),w=w[nmon][0].tolist())),flush=True)
 for beta in [0,.0001,.01,.1,.25,.5,.75,.9,.99,1.01,1.1,1.5,2,3,5,10,100,1000]:
  pp=np.exp(-beta*E);pp/=pp.sum(axis=1,keepdims=True)
  b=np.sum(w*(pp-p),axis=1);der=-np.sum(w*pp*(E-np.sum(pp*E,axis=1)[:,None]),axis=1);score=b*der*(1 if beta<1 else -1)
  good=(score>1e-8)&((vals[:,2]-vals[:,1])>1e-8)
  if np.any(good):
   z=np.argmax(np.where(good,score,-np.inf));out=dict(E=E[z].tolist(),rates=rates[z].tolist(),w=w[z].tolist(),vals=vals[z].tolist(),beta=beta,b=float(b[z]),der=float(der[z]),score=float(score[z]),count=count,accepted=accepted,seed=742908)
   print('CANDIDATE',json.dumps(out),flush=True);json.dump(out,open('research/adversarial_signed_candidate.json','w'),indent=2);raise SystemExit
 if batch%100==0:print(batch,count,accepted,nonmon,time.time()-start,flush=True)
json.dump(dict(count=count,accepted=accepted,nonmonotone_eigenvectors=int(nonmon),seed=742908),open('research/adversarial_signed_result.json','w'))
