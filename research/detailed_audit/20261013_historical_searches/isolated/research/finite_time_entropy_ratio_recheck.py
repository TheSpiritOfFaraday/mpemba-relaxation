"""80-digit check of the general-G contrast ratio-order obstruction."""
import json,mpmath as mp
from pathlib import Path
mp.mp.dps=80
c=json.load(open('research/finite_time_entropy_independent_audit.json'))['general_G']['first_ratio_order_failure']
n=c['n'];q=mp.matrix([mp.mpf(str(x)) for x in c['q']]);q/=sum(q)
M=mp.matrix(n)
for i in range(n):
 for j in range(i+1,n):
  k=mp.mpf(str(c['M'][i][j]))/q[i]
  M[i,j]=q[i]*k;M[j,i]=q[j]*k
for j in range(n):M[j,j]=-sum(M[i,j] for i in range(n) if i!=j)
E=mp.matrix([-mp.log(q[i]/q[0]) for i in range(n)]);b=mp.mpf(str(c['beta']));t=mp.mpf(str(c['physical_time']))
p=mp.matrix([mp.exp(-b*x) for x in E]);p/=sum(p);mu=sum(p[i]*E[i] for i in range(n))
r=mp.matrix([p[i]/q[i] for i in range(n)]);dr=mp.matrix([r[i]*(mu-E[i]) for i in range(n)])
P=mp.expm(M.T*t);rt=P*r;dt=P*dr
Z=mp.matrix(n,n-1)
for j in range(n-1):
 Q=sum(q[:j+1]);Qn=Q+q[j+1]
 for i in range(j+1):Z[i,j]=mp.sqrt(q[i]*q[j+1]/(Q*Qn))
 Z[j+1,j]=-mp.sqrt(Q/Qn)
x=Z.T*mp.matrix([mp.sqrt(q[i])*(rt[i]-1) for i in range(n)])
y=Z.T*mp.matrix([mp.sqrt(q[i])*dt[i] for i in range(n)])
ratios=[y[i]/x[i] for i in range(n-1)]
assert ratios[-1]>ratios[-2]
kl=sum(q[i]*dt[i]*mp.log(rt[i]) for i in range(n))
assert kl>0
out=dict(dps=80,ratios=[str(x) for x in ratios],last_adjacent_increase=str(ratios[-1]-ratios[-2]),KL_outward_derivative=str(kl),note='High precision numerical check, not an exact interval certificate; reconstructed reversible M from positive upper rates and normalized decimal q.')
Path('research/finite_time_entropy_ratio_recheck.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
