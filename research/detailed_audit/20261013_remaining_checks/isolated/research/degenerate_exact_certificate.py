"""Exact rank-two slow-space example and derivative-coefficient audit."""
import json,itertools
from pathlib import Path
import sympy as s
def require(ok,msg):
    if not ok:raise ValueError(msg)
x=s.Symbol('x');f=s.Matrix([1,x,x*x,x**3]);Z=sum(f)
q=s.Matrix([1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)]);p=q/sum(q);D=s.diag(*p)
ps=s.Matrix([1,s.Rational(3,4),s.Rational(9,16),s.Rational(27,64)]);ps/=sum(ps)
V=p.row_join(ps);H=D.inv()-D.inv()*V*(V.T*D.inv()*V).inv()*V.T*D.inv()
eps=s.Rational(1,100);M=-s.eye(4)+p*s.ones(1,4)+eps*D*H
require(H.T==H and H*D*H==H,'projector identities')
require(H.rank()==2 and H*p==s.zeros(4,1) and H*ps==s.zeros(4,1),'kernel and rank')
require(M*D==D*M.T and s.ones(1,4)*M==s.zeros(1,4),'reversibility and conservation')
require(all(M[i,j]>0 for i in range(4) for j in range(4) if i!=j),'positive rates')
l=s.Symbol('lambda');require(s.expand(M.charpoly(l).as_expr()-l*(l+1)*(l+s.Rational(99,100))**2)==0,'full spectrum')
Q=s.factor((f.T*H*f)[0]/Z**2);R=s.factor(-x*(s.diff((f.T*H*f)[0],x)*Z-2*(f.T*H*f)[0]*s.diff(Z,x)))
require(s.cancel(-x*s.diff(Q,x)-R/Z**3)==0,'derivative identity')
# Generic symbolic coefficients, independently expanded for arbitrary symmetric H.
ff=s.symbols('f1:5');ee=s.symbols('E1:5');HH=s.zeros(4)
for i in range(4):
    for j in range(i,4):HH[i,j]=HH[j,i]=s.Symbol('h%d%d'%(i,j))
W=(s.Matrix(ff).T*HH*s.Matrix(ff))[0]
RR=s.Poly(s.expand(-sum(ee[i]*ff[i]*s.diff(W,ff[i]) for i in range(4))*sum(ff)+2*W*sum(ee[i]*ff[i] for i in range(4))),*ff)
require(len(RR.terms())==16,'sixteen generic terms')
for i,k in itertools.permutations(range(4),2):
    mono=[0]*4;mono[i]=2;mono[k]=1
    require(s.expand(RR.coeff_monomial(tuple(mono))-2*(ee[k]-ee[i])*(HH[i,i]-HH[i,k]))==0,'double index coefficient')
for i,j,k in itertools.combinations(range(4),3):
    mono=[0]*4
    for v in (i,j,k):mono[v]=1
    target=2*(HH[i,j]*(2*ee[k]-ee[i]-ee[j])+HH[i,k]*(2*ee[j]-ee[i]-ee[k])+HH[j,k]*(2*ee[i]-ee[j]-ee[k]))
    require(s.expand(RR.coeff_monomial(tuple(mono))-target)==0,'distinct index coefficient')
numer=s.fraction(Q)[0]
require(numer.subs(x,s.Rational(1,2))==0 and numer.subs(x,s.Rational(3,4))==0,'two thermal zeros')
roots=s.Poly(numer,x).intervals(eps=s.Rational(1,10**12))
positive_roots=[str(r) for r in roots if r[0][1]>0]
out=dict(verified=True,scope='rank-two rational realization and generic sixteen-term derivative identity',
    bath_x='1/2',strong_x='3/4',spectrum=['0','-99/100','-99/100','-1'],
    generator=[[str(M[i,j]) for j in range(4)] for i in range(4)],Q=str(Q),R=str(R),
    positive_numerator_roots=positive_roots,generic_terms=len(RR.terms()),triplet_constraints_claimed=False)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
