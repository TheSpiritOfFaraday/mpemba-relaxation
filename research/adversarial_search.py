import numpy as np, json, time
from itertools import combinations
rng=np.random.default_rng(934256)
pairs=list(combinations(range(4),2)); trips=list(combinations(range(4),3))
start=time.time(); count=accepted=0
for batch in range(2000):
 n=10000
 E=np.column_stack([np.zeros(n),np.cumsum(10**rng.uniform(-3,1,(n,3)),axis=1)])
 q=np.exp(-E); p=q/q.sum(axis=1,keepdims=True)
 rates=10**rng.uniform(-5,5,(n,6)); rates/=rates.max(axis=1,keepdims=True)
 keep=np.ones(n,bool)
 for i,j,k in trips:
  A=rates[:,pairs.index((i,k))];B=rates[:,pairs.index((j,k))];C=rates[:,pairs.index((i,j))]
  t=q[:,j]/q[:,i];u=q[:,k]/q[:,j];h=A+B-C*(1+t)
  F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
  keep &= ~(((h>=0)&(F>0))|((h<=0)&(G>0)))
 E=E[keep];q=q[keep];p=p[keep];rates=rates[keep];m=len(E);count+=n;accepted+=m
 if not m:continue
 S=np.zeros((m,4,4));diag=np.zeros((m,4))
 for z,(i,j) in enumerate(pairs):
  a=rates[:,z];up=a*q[:,j]/q[:,i]
  S[:,i,j]=S[:,j,i]=a*np.sqrt(q[:,j]/q[:,i]);diag[:,j]-=a;diag[:,i]-=up
 for i in range(4): S[:,i,i]=diag[:,i]
 vals,vec=np.linalg.eigh(S);w=vec[:,:,2]/np.sqrt(p)
 # stable significance excludes extreme eigengaps and unstable numerical contamination
 for beta in [0,0.0001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1.01,1.1,1.5,2,3,5,10,100,1000]:
  pp=np.exp(-beta*E);pp/=pp.sum(axis=1,keepdims=True)
  b=np.sum(w*(pp-p),axis=1);der=-np.sum(w*pp*(E-np.sum(pp*E,axis=1)[:,None]),axis=1)
  score=b*der*(1 if beta<1 else -1)
  good=(score>1e-9)&((vals[:,2]-vals[:,1])>1e-10)&(vals[:,2]<-1e-10)
  if np.any(good):
   z=np.argmax(np.where(good,score,-np.inf));out=dict(E=E[z].tolist(),rates=rates[z].tolist(),w=w[z].tolist(),vals=vals[z].tolist(),beta=beta,b=float(b[z]),der=float(der[z]),score=float(score[z]),count=count,accepted=accepted,seed=934256)
   print('CANDIDATE',json.dumps(out),flush=True)
   json.dump(out,open('research/adversarial_candidate.json','w'),indent=2);raise SystemExit
 if batch%20==0: print(batch,count,accepted,'secs',time.time()-start,flush=True)
json.dump(dict(count=count,accepted=accepted,seed=934256),open('research/adversarial_result.json','w'))
