"""Exact rational certificate. All acceptance checks use SymPy rationals/integers.
Floats only locate intervals; correctness follows solely from exact tests.
"""
import sympy as s
from itertools import combinations
from fractions import Fraction
import json
q=list(map(s.Rational,['1','199/200','2/5','1/400']))
a=list(map(s.Rational,['1','28','249/10','1799/8','4983/200','251/25']))
pairs=list(combinations(range(4),2));M=s.zeros(4)
for v,(i,j) in zip(a,pairs): M[i,j]=v;M[j,i]=v*q[j]/q[i]
for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
assert all(sum(M[:,j])==0 for j in range(4))
assert all(M[i,j]>0 and M[i,j]*q[j]==M[j,i]*q[i] for i in range(4) for j in range(4) if i!=j)
x=s.Symbol('x');char=s.factor(M.charpoly(x).as_expr());poly=s.Poly(char/x,x);primitive=poly.clear_denoms()[1].primitive()[1]
print('GENERATOR',M);print('CHAR',char);print('PRIMITIVE_CUBIC',primitive.as_expr())
trip=[]
for i,j,k in combinations(range(4),3):
 A=M[i,k];B=M[j,k];C=M[i,j];t=q[j]/q[i];u=q[k]/q[j]
 H=A+B-C*(1+t);F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
 assert F<0 and G<0
 kappa=(A+B)/(C*(1+t))
 assert not ((kappa>=1 and F>0) or (kappa<=1 and G>0))
 trip.append(dict(triplet=[i+1,j+1,k+1],kappa=str((A+B)/(C*(1+t))),F=str(F),G=str(G)))
# Reset family alpha>24: M(alpha)=M(25)+(alpha-25)*(q*1^T-Z I).
# Each F and G is exactly invariant under this addition.
zreset=sum(q)
R=s.Matrix(q)*s.ones(1,4)-zreset*s.eye(4)
assert all(sum(R[:,j])==0 for j in range(4))
for i,j,k in combinations(range(4),3):
 t=q[j]/q[i];u=q[k]/q[j]
 assert s.simplify(q[i]*(1+t*u/2)-q[j]*(1+u/2)-q[i]*(1-t))==0
 assert q[i]*t-q[j]==0
print('TRIPLETS',trip)
# Exact isolating intervals produced by real-root isolation.
iso=s.polys.polytools.intervals(primitive,eps=s.Rational(1,10**25));print('ISOLATING_INTERVALS',iso)
assert len(iso)==3 and all(mult==1 and hi<0 for ((lo,hi),mult) in iso)
for ((lo,hi),mult) in iso: assert primitive.eval(lo)*primitive.eval(hi)<0
for z in range(2):assert iso[z][0][1]<iso[z+1][0][0]
L,U=iso[-1][0]
B=M.T-x*s.eye(4);co=[s.factor((-1)**j*B.minor_submatrix(0,j).det()) for j in range(4)]
assert all(s.rem(z,primitive.as_expr(),x)==0 for z in B*s.Matrix(co))
assert s.rem(sum(q[i]*co[i] for i in range(4)),primitive.as_expr(),x)==0
assert all(s.rem(z,primitive.as_expr(),x)==0 for z in (M-x*s.eye(4))*s.diag(*q)*s.Matrix(co))
print('COFACTORS',co)
# Tiny exact interval package. Outward correctness follows from rational arithmetic.
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def mul(a,b):
 t=[v*w for v in a for w in b];return(min(t),max(t))
def inv(a):
 assert a[0]*a[1]>0
 return(1/a[1],1/a[0])
def div(a,b):return mul(a,inv(b))
def peval(expr,I):
 out=(s.Rational(0),s.Rational(0))
 for c in s.Poly(expr,x).all_coeffs():out=add(mul(out,I),(c,c))
 return out
v=[peval(c,(L,U)) for c in co];w=[div(z,v[0]) for z in v]
print('LEFT_EIGENVECTOR_INTERVALS',[(str(l),str(u)) for l,u in w])
# v = diag(q) w / sum_i q_i w_i**2 is the right mode and w.T v = 1.
# The preceding exact right residual proves its eigenvalue equation.
normw=(s.Rational(0),s.Rational(0))
for qi,wi in zip(q,w):normw=add(normw,mul((qi,qi),mul(wi,wi)))
assert normw[0]>0
# Certify q**(1/n) by integer/rational powers. Locate via decimal approximation,
# but root enclosure acceptance and final inequalities use exact arithmetic.
def root_bound(r,n,digits=24):
 den=10**digits
 guess=int(s.floor(s.N(r**s.Rational(1,n),digits+15)*den))
 lo=s.Rational(guess,den);hi=s.Rational(guess+1,den)
 assert lo>=0 and lo**n<=r<=hi**n
 return lo,hi

def overlap(n):
 z=[root_bound(r,n) for r in q];norm=(s.Rational(0),s.Rational(0));num=(s.Rational(0),s.Rational(0))
 for zi,wi in zip(z,w):norm=add(norm,zi);num=add(num,mul(zi,wi))
 return div(num,norm),z
hot,zh=overlap(1000);cold,zc=overlap(100)
assert cold[1]<hot[0]<hot[1]<0
margin=add(hot,neg(cold));assert margin[0]>0
print('HOT_B',tuple(s.N(z,30) for z in hot));print('COLD_B',tuple(s.N(z,30) for z in cold));print('MARGIN',tuple(s.N(z,30) for z in margin));print('ALL EXACT CERTIFICATE CHECKS PASSED')
report=dict(q=list(map(str,q)),down_rates=list(map(str,a)),generator=[[str(M[i,j]) for j in range(4)] for i in range(4)],cubic=str(primitive.as_expr()),eigenvalue_intervals=[list(map(str,I)) for I,m in iso],cofactor_polynomials=list(map(str,co)),left_intervals=[[str(z) for z in I] for I in w],triplets=trip,beta_hot='1/1000',beta_cold='1/100',hot_root_intervals=[[str(z) for z in I] for I in zh],cold_root_intervals=[[str(z) for z in I] for I in zc],hot_overlap_interval=[str(z) for z in hot],cold_overlap_interval=[str(z) for z in cold],overlap_margin_lower_bound=str(margin[0]),checks='passed')
json.dump(report,open('research/adversarial_certificate.json','w'),indent=2)
