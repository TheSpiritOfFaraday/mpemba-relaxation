"""Exact identities for a physically reachable failure of a proof-route ratio order.
No claim of a KL reversal is made; finite beta/time existence is analytic.
"""
import sympy as s,itertools,json
from pathlib import Path
p,q,r,h,a=s.symbols('p q r h a',positive=True);weights=s.Matrix([p,q,r,h]);D=s.diag(*weights)
pairs=list(itertools.combinations(range(4),2));base=[2,2,1,2,2,3];M=s.zeros(4);checks=[]
for (i,j),v in zip(pairs,base):M[i,j]=weights[i]*v+a;M[j,i]=weights[j]*(v+a/weights[i])
for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
expectedF=[0,-p-h/2,-p-r-h,-r-h/2];expectedG=[0,-q,-2*r,-r]
for (i,j,k),F0,G0 in zip(itertools.combinations(range(4),3),expectedF,expectedG):
 A,B,C=M[i,k],M[j,k],M[i,j];t=weights[j]/weights[i];u=weights[k]/weights[j]
 F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
 assert s.factor(F-F0-a*(t-1)*(1+u/2))==0
 assert s.factor(G-G0-a*(t-1))==0
 checks.extend([f'F{i+1}{j+1}{k+1} strict-repair identity',f'G{i+1}{j+1}{k+1} strict-repair identity'])
H=s.Matrix([[1,1,1],[-p/q,1,1],[0,-(p+q)/r,1],[0,0,-(p+q+r)/h]])
T=H.T*D*M.T*H
assert s.factor(T[2,0]-T[1,0]-p*(p+q+r))==0
checks.append('sqrt(eta1) times f-prime(0)=q1 Q3>0')
out={'verified':True,'checks':checks,'scope':'Exact all-strict triplet identities and positive initial-time derivative of the cold-preparation limit. The companion analytic argument establishes finite thermal reachability; no KL reversal claim.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
