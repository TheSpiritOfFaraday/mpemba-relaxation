"""Independent state-space checks of the weighted-contrast theorem.
Numerical stress checks supplement, rather than prove, the all-n theorem.
"""
import json,time,itertools
from pathlib import Path
import numpy as np
from scipy.linalg import expm
rng=np.random.default_rng(111031);start=time.time();counts=0;errors=[];min_order=1.;max_identity=0.
for n in [3,4,5,7,10]:
 for case in range(200):
    E=np.r_[0,np.cumsum(rng.uniform(.05,1.5,n-1))];q=np.exp(-E);q/=sum(q);D=np.diag(q)
    K=np.zeros((n,n))
    for j in range(1,n):
        K[:j,j]=np.sort(np.exp(rng.uniform(-4,4,j)));K[j,:j]=K[:j,j]
    L=K*q[None,:];L-=np.diag(L.sum(axis=1));M=L.T
    h=np.zeros((n,n-1))
    for k in range(n-1):
        Q=sum(q[:k+1]);eta=Q*sum(q[:k+2])/q[k+1]
        h[:k+1,k]=1/np.sqrt(eta);h[k+1,k]=-Q/q[k+1]/np.sqrt(eta)
    A=h.T@D@L@h;off=A-np.diag(np.diag(A))
    if off.min()<-1e-9:errors.append(['negative contrast entry',n,case,float(off.min())])
    identity=max(np.max(abs(h.T@D@h-np.eye(n-1))),np.max(abs(L@h-h@A)))
    max_identity=max(max_identity,float(identity))
    # Reconstruct the original state-space slow projector independently.
    S=M*np.sqrt(q[None,:]/q[:,None]);vals,V=np.linalg.eigh(S);v=V[:,-2]
    def p(beta):
        z=np.exp(-beta*E);return z/sum(z)
    side_sign=[]
    for beta in [.01,.3,.8,1.2,2.,10.]:
        b=v@(p(beta)/np.sqrt(q));side_sign.append(b)
    if not(all(x*y>0 for x,y in zip(side_sign[:2],side_sign[1:3])) and all(x*y>0 for x,y in zip(side_sign[3:5],side_sign[4:])) and side_sign[2]*side_sign[3]<0):errors.append(['slow overlap sign',n,case,side_sign])
    # Original p(t)=exp(tM)p, rather than a contrast-coordinate simulation.
    near,far=p(1.2),p(2.5)
    for scaled_time in [0,.01,.1,1.,3.]:
        t=scaled_time/(-vals[-2]);U=expm(M*t);pn=U@near;pf=U@far
        dn=sum((pn-q)**2/q);df=sum((pf-q)**2/q)
        margin=(df-dn)/max(df,dn,1e-300);min_order=min(min_order,float(margin));counts+=1
        if not df>dn:errors.append(['heating chi-square order',n,case,t,dn,df])
out=dict(seed=111031,systems=1000,state_space_heating_checks=counts,errors=errors,
    minimum_relative_heating_margin=min_order,max_contrast_identity_residual=max_identity,
    elapsed_seconds=time.time()-start,scope='Floating-point original-state-space verification in dimensions 3,4,5,7,10')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if errors:raise SystemExit('Audit failures require investigation')
