"""Exact algebra audit for the sharp four-state full-triplet amplitude bound.

Checks the Helmert-to-Schur identities and all three normalized-face
positivity decompositions. Symbolic equalities support the independent
analytic sign, spectral and limiting arguments in the TeX proof.
"""
from pathlib import Path
import json,time
import sympy as s
start=time.time();checks=[]
def require(condition,label):
    if not condition:raise AssertionError(label)
    checks.append(label)
def zero(expr,label):require(s.factor(expr)==0,label)
p,q,r,h=s.symbols('p q r h',positive=True)
a,b,c,d,e,f=s.symbols('a b c d e f',nonnegative=True)
weights=[p,q,r,h];rates=[a,b,c,d,e,f];pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
Q2=p+q;Q3=p+q+r;Q4=p+q+r+h
M=s.zeros(4)
for (i,j),rate in zip(pairs,rates):
    M[i,j]=rate;M[j,i]=rate*weights[j]/weights[i]
for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
H=s.Matrix([[1,1,1],[-p/q,1,1],[0,-Q2/r,1],[0,0,-Q3/h]])
D=s.diag(*weights);eta=s.diag(p*Q2/q,Q2*Q3/r,Q3*Q4/h)
require(s.simplify(H.T*D*H-eta)==s.zeros(3),'weighted Helmert normalization')
B=s.simplify(eta.inv()*H.T*D*M.T*H)
F123=b*(1+r/(2*p))-d*(1+r/(2*q))-a*(1-q/p)
F124=c*(1+h/(2*p))-e*(1+h/(2*q))-a*(1-q/p)
G123=b*q/p-d;G124=c*q/p-e;G134=c*r/p-f;G234=e*r/q-f
ss=(p-q)*r/(2*q*Q3)
C=s.factor(B[1,0]+ss*(B[1,1]-B[0,0])-ss**2*B[0,1])
DS=s.factor((B[1,2]-ss*B[0,2])*(B[2,0]+ss*B[2,1]))
zero(B[1,2]-ss*B[0,2]+Q4*(G134+G234)/(2*Q3),'first Schur cross factor is nonnegative by G')
N=-2*p*q*r*F123;J=e*p*(2*q+r)-c*q*(2*p+r)
T=(c+e+f)*N+h*(2*c*e*r*(p-q)+f*J)
pref=Q2*Q4/(4*p*q*q*Q3*Q3)
zero(-B[2,2]*C+DS-pref*T,'eigenvalue-free Schur expression factorization')
lam=s.Symbol('lambda');delta=lam-B[2,2]
E=s.Matrix(2,2,lambda i,j:B[i,j]+B[i,2]*B[2,j]/delta)
zero(delta*(E[1,0]+ss*(E[1,1]-E[0,0])-ss**2*E[0,1])-delta*C-DS,'Schur quadratic at candidate ratio')
# Face c=0: all remaining displayed summands are nonnegative.
zero(T.subs(c,0)-((e+f)*N+h*f*e*p*(2*q+r)),'face a14=0 positivity decomposition')
# Face a=0: J is a positive multiple of -F124 plus a positive remainder.
Jface=(2*p+r)/(2*p+h)*(-2*p*q*F124.subs(a,0))+2*e*p*(p-q)*(r-h)/(2*p+h)
zero(J-Jface,'face a12=0 positivity decomposition')
# Face b=0: N is manifestly nonnegative, and the f coefficient is controlled.
Nface=2*a*q*r*(p-q)+d*p*r*(2*q+r)
zero(N.subs(b,0)-Nface,'face a13=0 numerator identity')
coef=Nface+h*J
controlled=2*p*q*(r-h)*(c-e)+d*p*r*(2*q+r)-2*p*q*r*F124
zero(coef-controlled,'face a13=0 coefficient controlled by F124')
zero(T.subs(b,0)-((c+e)*Nface+f*coef+2*c*e*r*h*(p-q)),'face a13=0 total decomposition')
zero(J-(e-c)*p*(2*q+r)-c*r*(p-q),'J nonnegative whenever c<=e')
# Reset translation preserves every linear test and shifts B by a scalar.
tau=s.Symbol('tau');reset={a:a-tau*p,b:b-tau*p,c:c-tau*p,d:d-tau*q,e:e-tau*q,f:f-tau*r}
for name,expr in [('F123',F123),('F124',F124),('G123',G123),('G124',G124),('G134',G134),('G234',G234)]:zero(expr.subs(reset,simultaneous=True)-expr,'reset invariance '+name)
require(s.simplify(B.subs(reset,simultaneous=True)-B-tau*Q4*s.eye(3))==s.zeros(3),'reset shift of contrast generator')
# Conditional three-state monotonicity identity.
x,y=s.symbols('x y',positive=True);aa,bb=s.symbols('aa bb',positive=True)
Z=1+x+y;fun=(1-x)/Z
calc=s.diff(fun,x)*(-aa*x)+s.diff(fun,y)*(-bb*y)
zero(calc-(2*aa*x+bb*y-(bb-aa)*x*y)/Z**2,'conditional first contrast derivative')
# The numerator is >0 for 0<x<=1 and 0<aa<bb.
zero(2*aa*x+bb*y-(bb-aa)*x*y-(2*aa*x+aa*x*y+bb*y*(1-x)),'conditional derivative positive decomposition')

out={'verified':True,'exact_checks':checks,'elapsed_seconds':time.time()-start,
 'assumptions':'p>q>r>h>0; shifted rates nonnegative; all G<=0; F123,F124<=0; at least one of a,b,c=0',
 'scope':'Exact symbolic identities and manifest-sign decompositions; spectral Schur arguments and sharpness are proved in full_triplet_amplitude_bound.tex',
 'additional_F_needed':['F123','F124'],'amplitude_ratio_bound':'strictly less than 4/3','amplitude_advantage_supremum':'1/4'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
