import sympy as s
from itertools import combinations
p,q,r,h=s.symbols('p q r h',positive=True);bb,cc=s.symbols('b c',positive=True)
qs=s.Matrix([p,q,r,h]);D=s.diag(*qs);Q=p+q+r+h
H=s.Matrix([[1,1,1],[-p/q,1,1],[0,-(p+q)/r,1],[0,0,-(p+q+r)/h]])
N=s.diag(p*(p+q)/q,(p+q)*(p+q+r)/r,(p+q+r)*Q/h);v=s.Matrix([1,bb,cc]);ell=s.symbols('l m n',nonnegative=True);K=s.zeros(3)
for (i,j),e in zip(combinations(range(3),2),ell):K[i,j]=e/N[i,i];K[j,i]=e/N[j,j]
for i in range(3):K[i,i]=-sum(K[i,j]*v[j] for j in range(3) if i!=j)/v[i]
M=(H*K*N.inv()*H.T*D).T
u,v,m=s.symbols('u v m',nonnegative=True)
for i,j,k in [(0,1,2),(0,1,3)]:
 t=qs[j]/qs[i];uu=qs[k]/qs[j]
 FF=M[i,k]*(1+t*uu/2)-M[j,k]*(1+uu/2)-M[i,j]*(1-t)
 FF=FF.subs(dict(zip(ell,[u+h/Q*m,m,m+v])), simultaneous=True)
 print('F',i,j,k,[s.factor(s.diff(FF,z)) for z in [u,v,m]],flush=True)
