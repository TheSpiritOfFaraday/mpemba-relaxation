"""Exploratory symbolic realization; numeric diagnostics are not certification."""
import sympy as s
from itertools import combinations
z=s.symbols('z',real=True); q=s.Matrix([1,s.Rational(199,200),s.Rational(2,5),s.Rational(1,400)])
w=s.Matrix([1,s.Rational(-759,1000),0,z]); w[2]=s.factor(-sum(q[i]*w[i] for i in [0,1,3])/q[2])
x,y=s.symbols('x y'); a=[s.Rational(1),s.Rational(28),s.Rational(249,10),x,s.Rational(4983,200),y]; M=s.zeros(4)
for (i,j),v in zip(combinations(range(4),2),a): M[i,j]=v;M[j,i]=v*q[j]/q[i]
for j in range(4): M[j,j]=-sum(M[:,j])
lam=s.factor((M.T*w)[0]); sol={x:s.factor((lam*w[1]-a[0]*(1-w[1])-a[4]*q[3]/q[1]*(z-w[1]))/(q[2]/q[1]*(w[2]-w[1]))),y:s.factor((lam*z-a[2]*(1-z)-a[4]*(w[1]-z))/(w[2]-z))}; M=M.subs(sol); a=[v.subs(sol) for v in a]
print('w',w.T,'lambda',lam,'rates',a)
for i,j,k in combinations(range(4),3):
 A=M[i,k];B=M[j,k];C=M[i,j];t=q[j]/q[i];u=q[k]/q[j]
 F=s.factor(A*(1+t*u/2)-B*(1+u/2)-C*(1-t));G=s.factor(A*t-B)
 print('trip',i+1,j+1,k+1,'F=',F,'G=',G)
E=s.Matrix([-s.log(qi) for qi in q]); N0=sum(w); cov=s.expand(sum(w[i]*E[i] for i in range(4))-sum(w)*sum(E)/4); root=-cov.subs(z,0)/s.diff(cov,z)
print('N0',N0,'4Cov',cov,'zcrit',root,s.N(root,20))
for zz in [-.1,-.02,-.01,-.008,-.006,-.004,0,.005,.01,.02,.1]:
 mat=M.subs(z,zz);import numpy as np
 print(zz,'rates',[round(float(v.subs(z,zz)),7) for v in a],'eigs',np.linalg.eigvals(np.array(mat,dtype=float)))
