"""Exact counterexample to repeated-slow => stochastic-monotone in n=5."""
import sympy as s,itertools,json
from pathlib import Path
q=s.Matrix([s.Rational(v,43) for v in [24,12,4,2,1]])
r=s.Matrix([s.Rational(3,5),s.Rational(3,10),s.Rational(1,10)]);z=s.Matrix([1,2,-12]);D3=s.diag(*r)
J3=r*s.ones(1,3)-s.eye(3)-D3*z*z.T/(z.T*D3*z)[0]
J=s.zeros(5);J[:3,:3]=J3;M=q*s.ones(1,5)-s.eye(5)+J
assert all(M[i,j]>0 for i in range(5) for j in range(5) if i!=j)
assert s.ones(1,5)*M==s.zeros(1,5)
assert M*s.diag(*q)==s.diag(*q)*M.T
x=s.Symbol('x');assert s.factor(M.charpoly(x).as_expr())==x*(x+1)**2*(x+2)*(x+3)
table=[]
for i,j,k in itertools.combinations(range(5),3):
 A,B,C=M[i,k],M[j,k],M[i,j]
 F=s.factor(A*(1+q[k]/(2*q[i]))-B*(1+q[k]/(2*q[j]))-C*(1-q[j]/q[i]));G=s.factor(A*q[j]/q[i]-B)
 assert F<=0 and G<=0
 table.append({'triplet':[i+1,j+1,k+1],'F':str(F),'G':str(G)})
K21=s.factor(M[0,1]-M[0,2]);assert K21==-s.Rational(14,27)
result={'verified':True,'q':[str(v) for v in q],'generator':[[str(v) for v in row] for row in M.tolist()],'spectrum':['0','-1','-1','-2','-3'],'K21':str(K21),'triplets':table,'scope':'Repeated slow eigenspace does not imply stochastic monotonicity in five states. No claim of a Mpemba witness.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
