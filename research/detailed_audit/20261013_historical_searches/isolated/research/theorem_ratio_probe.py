"""Adversarial numerical audit of a proposed Schur-complement certificate, not a proof."""
import numpy as np,sys,json
from cone_search import rays,generator
rng=np.random.default_rng(115103);pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
negcount=0;minratio=1e9
for case in range(300):
 E=np.r_[0,np.cumsum(np.exp(rng.uniform(-3,2,3)))];q=np.exp(-E);Q=q.sum();rs=rays(q)
 if len(rs)==0:continue
 rates=np.exp(rng.uniform(-8,4,(300,len(rs))))@rs
 reset=np.array([q[i] for i,j in pairs]);rates-=np.min(rates/reset,axis=1)[:,None]*reset
 H=np.array([[1,1,1],[-q[0]/q[1],1,1],[0,-(q[0]+q[1])/q[2],1],[0,0,-sum(q[:3])/q[3]]]);etas=np.diag(H.T@np.diag(q)@H)
 B=np.array([(H.T@np.diag(q)@generator(q,a).T@H)/etas[:,None] for a in rates])
 ss=(q[0]-q[1])*q[2]/(2*q[1]*sum(q[:3]));C=B[:,1,0]+ss*(B[:,1,1]-B[:,0,0])-ss**2*B[:,0,1];D=(B[:,1,2]-ss*B[:,0,2])*(B[:,2,0]+ss*B[:,2,1]);T=-B[:,2,2]*C+D
 negcount+=int(np.sum(C < -1e-8));minratio=min(minratio,float(np.min(T/(np.max(abs(B),axis=(1,2))**2+1e-99))))
 bad=np.flatnonzero((C < -1e-8)*(T < -1e-7))
 if len(bad):
  j=bad[0];print(json.dumps({'case':case,'q':q.tolist(),'rates':rates[j].tolist(),'C':C[j],'T_bound':T[j],'B':B[j].tolist(),'s':ss}),flush=True);break
else:print('none',negcount,minratio)
