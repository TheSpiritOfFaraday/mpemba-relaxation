"""Construct near-boundary rational input numerically, decide it exactly."""
import json,time
from pathlib import Path
import sympy as s
import mpmath as mp
from exact_rational_triplet_classifier import classify,PAIRS
mp.mp.dps=80;start=time.time()
q0=list(map(s.Rational,['1','199/200','2/5','1/400']));q=[z/sum(q0) for z in q0]
r0=list(map(s.Rational,['1','28','249/10','1799/8','4983/200','251/25']))
conv=lambda z:mp.mpf(int(z.p))/int(z.q)
qm=list(map(conv,q));E=[-mp.log(z) for z in qm]
def numeric_endpoint(t):
 rates=[1+t*(z-1) for z in r0];M=mp.matrix(4)
 for (i,j),r in zip(PAIRS,rates):M[i,j]=conv(r);M[j,i]=conv(r*q[j]/q[i])
 for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
 S=mp.matrix([[M[i,j]*mp.sqrt(qm[j]/qm[i]) for j in range(4)] for i in range(4)])
 eig,U=mp.eigsy(S);w=[U[i,2]/mp.sqrt(qm[i]) for i in range(4)]
 return sum(w)/4*(-sum(w[i]*(E[i]-sum(E)/4) for i in range(4))/4)
lo=s.Rational(0);hi=s.Rational(1)
assert numeric_endpoint(lo)<0<numeric_endpoint(hi)
for _ in range(100):
 mid=(lo+hi)/2
 if numeric_endpoint(mid)>0:hi=mid
 else:lo=mid
results=[]
for label,t in [('below',lo),('above',hi)]:
 rates=[1+t*(z-1) for z in r0]
 bounded=classify(q0,rates,max_rounds=1)
 assert bounded['verified'] is False and bounded['ordinary_cooling'] is None and bounded['endpoint_sign'] is None
 exact=classify(q0,rates)
 assert exact['verified'] and exact['ordinary_cooling']==(label=='above')
 results.append(dict(label=label,t=str(t),weights=list(map(str,q0)),rates=list(map(str,rates)),bounded=bounded,exact=exact))
out=dict(verified=True,elapsed_seconds=time.time()-start,note='Floating eigensolver used only to construct candidates; all accepted signs and inconclusive outcomes are exact classifier outputs.',results=results)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(verified=True,elapsed_seconds=out['elapsed_seconds'],cases=[dict(label=z['label'],bounded_verified=z['bounded']['verified'],exact_cooling=z['exact']['ordinary_cooling'],exact_rounds=z['exact']['interval_rounds']) for z in results]),indent=2))
