import sympy as s
from itertools import combinations
p,q,r,h=s.symbols('p q r h',positive=True);bb,cc=s.symbols('b c',positive=True)
qs=s.Matrix([p,q,r,h]);D=s.diag(*qs);Q=p+q+r+h
H=s.Matrix([[1,1,1],[-p/q,1,1],[0,-(p+q)/r,1],[0,0,-(p+q+r)/h]])
N=s.diag(p*(p+q)/q,(p+q)*(p+q+r)/r,(p+q+r)*Q/h);v=s.Matrix([1,bb,cc]);ell=s.symbols('l m n',nonnegative=True);K=s.zeros(3)
for (i,j),e in zip(combinations(range(3),2),ell):K[i,j]=e/N[i,i];K[j,i]=e/N[j,j]
for i in range(3):K[i,i]=-sum(K[i,j]*v[j] for j in range(3) if i!=j)/v[i]
M=(H*K*N.inv()*H.T*D).T
for i,j,k in combinations(range(4),3):
 t=qs[j]/qs[i];u=qs[k]/qs[j]
 F=s.factor(M[i,k]*(1+t*u/2)-M[j,k]*(1+u/2)-M[i,j]*(1-t));G=s.factor(M[i,k]*t-M[j,k]);print('TRIP',i,j,k)
 for expr in [F,G]:print([s.factor(s.diff(expr,z)) for z in ell])
