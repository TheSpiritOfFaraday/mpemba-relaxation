import sympy as s
from itertools import combinations
q=s.symbols('p q r h',positive=True);a=s.symbols('a b c d e f');pairs=list(combinations(range(4),2));L=s.zeros(4)
for (i,j),v in zip(pairs,a):L[j,i]=v;L[i,j]=v*q[j]/q[i]
for i in range(4):L[i,i]=-sum(L[i,:])
T=s.Matrix([[1,-1,0,0],[1,1,-2,0],[q[0],q[1],q[2],-sum(q[:3])],list(q)])
K=s.simplify(T*L*T.inv())
for i in range(3):
 for j in range(3):
  if i!=j:print(i,j,s.factor(K[i,j]))
