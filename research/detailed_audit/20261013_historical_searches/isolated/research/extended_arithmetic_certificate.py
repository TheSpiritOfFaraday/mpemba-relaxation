import sympy as s
from itertools import combinations
p,q,r,h=s.symbols('p q r h',positive=True);bb,cc=s.symbols('b c',positive=True)
qs=s.Matrix([p,q,r,h]);D=s.diag(*qs);Q=p+q+r+h
H=s.Matrix([[1,1,1],[-p/q,1,1],[0,-(p+q)/r,1],[0,0,-(p+q+r)/h]])
N=s.diag(p*(p+q)/q,(p+q)*(p+q+r)/r,(p+q+r)*Q/h);v=s.Matrix([1,bb,cc]);ell=s.symbols('l m n',nonnegative=True);K=s.zeros(3)
for (i,j),e in zip(combinations(range(3),2),ell):K[i,j]=e/N[i,i];K[j,i]=e/N[j,j]
for i in range(3):K[i,i]=-sum(K[i,j]*v[j] for j in range(3) if i!=j)/v[i]
M=(H*K*N.inv()*H.T*D).T

# All comparisons below are symbolic identities. Positivity arguments are in
# extended_adversarial_progress.txt and use p>q>r>h>0 and 0<b<s.
u,v,m=s.symbols('u v m',nonnegative=True)
Q2=p+q;Q3=p+q+r;Q4=p+q+r+h
Fs=[];Gs=[]
for i,j,k in combinations(range(4),3):
 t=qs[j]/qs[i];uu=qs[k]/qs[j]
 Fs.append(s.factor(M[i,k]*(1+t*uu/2)-M[j,k]*(1+uu/2)-M[i,j]*(1-t)))
 Gs.append(s.factor(M[i,k]*t-M[j,k]))
expectedG=[q/(p*Q3)*(-ell[0]+h/Q4*ell[1]), -q/(p*Q4)*ell[1], -r/(Q2*Q4)*(q/p*ell[1]+ell[2]), r/(Q2*Q4)*(ell[1]-ell[2])]
assert all(s.factor(x-y)==0 for x,y in zip(Gs,expectedG))
subs=dict(zip(ell,[u+h/Q4*m,m,m+v]))
f=s.factor(Fs[0].subs(subs, simultaneous=True));g=s.factor(Fs[1].subs(subs, simultaneous=True))
A=-(bb*(p-q)+2*p)*(2*bb*q*Q3-(p-q)*r)/(2*bb*p*Q2**2*Q3)
B=-(p-q)*(bb*q*Q3-p*r)*(bb*h+cc*Q4+h)/(bb*p*Q2**2*Q3*Q4)
C=(p-q)*(cc*r*Q4-bb*h*Q2)/(bb*Q2**2*Q3*Q4)
assert s.factor(f-A*u-B*m-C*v)==0
D=-(p-q)*(bb*Q3+r)*(bb*q*Q3-p*r)/(bb*p*Q2**2*Q3**2)
k=(bb*Q2-cc*r)/(cc*Q3)
assert s.factor(s.diff(g,u)-D)==0
assert s.factor(s.diff(g,v)+k*C)==0
combo=s.factor(g+k*f)
assert s.factor(s.diff(combo,v))==0
assert s.factor(s.diff(combo,u)-(D+k*A))==0
Nm=s.factor(s.diff(combo,m)*2*cc*p*Q2*Q3*Q4)
assert s.factor(s.diff(Nm,cc,2)/2+2*q*(p-q)*Q4)==0
N0=2*h*(bb+1)*(p-q)*(p-bb*q)
assert s.factor(Nm.subs(cc,0)-N0)==0
cmax=bb*h*Q2/(r*Q4);Nend=s.factor(Nm.subs(cc,cmax));st=(p-q)*r/(2*q*Q3)
assert s.factor(s.diff(Nend,bb,2)/2+2*h*q*(h+r)*(p-q)*Q3**2/(r**2*Q4))==0
assert s.factor(Nend.subs(bb,0)-2*h*p*(p-q))==0
assert s.factor(Nend.subs(bb,st)-h*r*(p-q)*Q2**2/(2*q*Q4))==0
# Reconstruct both original-rate generators, conservation and reversibility.
assert s.simplify(s.ones(1,4)*M)==s.zeros(1,4)
assert s.simplify(M*s.diag(*qs)-s.diag(*qs)*M.T)==s.zeros(4)
assert s.simplify(M.T*H*s.Matrix([1,bb,cc]))==s.zeros(4,1)
print('ALL SYMBOLIC IDENTITIES PASSED')
print('A=',A);print('B=',B);print('C=',C);print('D=',D);print('k=',k)
print('Nm=',Nm)
