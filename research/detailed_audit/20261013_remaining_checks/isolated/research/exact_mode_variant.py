import sympy as s,numpy as np,sys
from cone_search import constraints,generator
q=s.Matrix([1,s.Rational(199,200),s.Rational(2,5),s.Rational(1,400)])
# Design a rational mode close to the independent signed-cone witness.
w=s.Matrix([1,s.Rational(-759,1000),0,s.Rational(-4,1000)])
w[2]=-(q[0]*w[0]+q[1]*w[1]+q[3]*w[3])/q[2]
a12,a13,a14,a23=map(s.Rational,['1','28','24.9','224.875'])
lam=(a12*q[1]*(w[1]-1)+a13*q[2]*(w[2]-1)+a14*q[3]*(w[3]-1))
x,y=s.symbols('x y');aa=[a12,a13,a14,x,s.Rational(4983,200),y];M=s.zeros(4)
from itertools import combinations
for (i,j),a in zip(combinations(range(4),2),aa):M[i,j]=a;M[j,i]=a*q[j]/q[i]
for j in range(4):M[j,j]=-sum(M[:,j])
sol=s.solve(list(M.T*w-lam*w),(x,y));M=M.subs(sol);aa=[a.subs(sol) for a in aa]
print('q',list(q),'w',list(w),'lambda',lam,'rates',aa)
print('constraints',constraints(np.array(q,dtype=float).ravel())@np.array(aa,dtype=float))
print('eigen',np.linalg.eigvals(np.array(M,dtype=float)))
print('charpoly',s.factor(M.charpoly().as_expr()))
for b in [.001,.01]:
 p=np.array(q,dtype=float).ravel()**b;p/=sum(p);print('b',b,float(np.array(w,dtype=float).ravel()@p))
