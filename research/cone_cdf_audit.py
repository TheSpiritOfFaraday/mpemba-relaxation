"""Numerical audit of slow-mode sign patterns inside the triplet-failure cone."""
import numpy as np,json,time
from pathlib import Path
from collections import Counter
from cone_search import rays,generator
rng=np.random.default_rng(111027);counter=Counter();found=[];start=time.time()
for case in range(500):
 E=np.r_[0,np.cumsum(10**rng.uniform(-2.5,.8,3))];q=np.exp(-E)
 rs=rays(q)
 if len(rs)==0:continue
 a=np.exp(rng.uniform(-8,3,(500,len(rs))))@rs
 reset=np.array([q[i] for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]);rho=np.maximum(0,(-a/reset).max(axis=1))+1;a+=rho[:,None]*reset
 M=np.array([generator(q,x) for x in a]);vals,V=np.linalg.eigh(M*np.sqrt(q[None,None,:]/q[None,:,None]));w=V[:,:,-2]/np.sqrt(q);w*=np.where(w[:,0]<0,-1,1)[:,None]
 for i,ww in enumerate(w):
  norm=abs(ww).max();sg=np.where(ww>1e-8*norm,1,np.where(ww< -1e-8*norm,-1,0));pattern=tuple(sg);counter[str(pattern)]+=1
  cum=np.cumsum(q*ww)[:3];counter['CDF_nonnegative' if min(cum)>=-1e-8*max(abs(cum)) else 'CDF_mixed']+=1
  changes=sum(sg[j]*sg[j+1]<0 for j in range(3))
  if min(cum)<-1e-8*max(abs(cum)) and len(found)<10:found.append(dict(E=E.tolist(),rates=a[i].tolist(),w=ww.tolist(),eigenvalues=vals[i].tolist(),pattern=list(map(int,sg))))
 if case%50==0:print('progress',case,sum(counter.values()),len(found),flush=True)
out={'seed':111027,'count':sum(v for k,v in counter.items() if k.startswith('(')),'patterns':dict(counter),'multiple_sign_change_examples':found,'elapsed':time.time()-start,'scope':'Numerical observations only; near-zero coordinates use relative1e-8 threshold'}
Path('research/cone_cdf_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
