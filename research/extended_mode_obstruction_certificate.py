"""Exact Farkas certificate refuting sufficiency of the arithmetic mode cone."""
import sympy as s,json
from itertools import combinations
def require(ok):
 if not ok:raise ValueError("exact obstruction identity failed")
q=s.Matrix([5,4,2,1]);w=s.Matrix([162,-18,71,-880]);D=s.diag(*q)
require((q.T*w)[0]==0 and w[0]>w[1] and w[0]+w[1]>2*w[2] and w[3]<0)
H=s.Matrix([[1,1,1],[-s.Rational(5,4),1,1],[0,-s.Rational(9,2),1],[0,0,-11]])
require(H*s.Matrix([80,2,80])==w)
U=s.Matrix.hstack(*(q.T.col_join((D*w).T)).nullspace())
b11,b12,b22=s.symbols('b11 b12 b22');B=s.Matrix([[b11,b12],[b12,b22]])
M=-D*U*B*U.T
require(s.ones(1,4)*M==s.zeros(1,4))
require(M.T*w==s.zeros(4,1))
vals=[]
for i,j,k in combinations(range(4),3):
 t=q[j]/q[i];u=q[k]/q[j]
 F=M[i,k]*(1+t*u/2)-M[j,k]*(1+u/2)-M[i,j]*(1-t)
 G=M[i,k]*t-M[j,k]
 vals.extend([s.factor(F),s.factor(G)])
weights=[s.Rational(482571629,11676414),s.Rational(496382581,11676414),s.Rational(16620682809,77842760)]
require(min(weights)>0)
require(s.factor(weights[0]*vals[0]+weights[1]*vals[2]+weights[2]*vals[7]-s.trace(B))==0)
# F123,F124,G234<=0 force trace(B)<=0. But any simple slow realization
# requires B positive definite by exact_mode_realizability theorem, impossible.
out=dict(q=list(map(str,q)),w=list(map(str,w)),helmert_coefficients=['80','2','80'],
         U=[[str(x) for x in U.row(i)] for i in range(4)],
         dual_constraints=['F123','F124','G234'],dual_weights=list(map(str,weights)),
         identity='positive weighted sum of F123,F124,G234 equals trace(B)',
         verified=True,conclusion='No full-triplet simple slow realization; strict necessary cone is not sufficient')
open('research/extended_mode_obstruction_certificate.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('EXACT MODE OBSTRUCTION CERTIFICATE PASSED')
