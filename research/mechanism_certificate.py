"""Exact certificate for a predictive one-parameter four-state family.

All acceptance tests use integer/rational arithmetic. Bernstein coefficients
certify polynomial signs on a CONTINUOUS interval, not a parameter grid.
Logarithms use an explicit convergent-series remainder; fractional Gibbs powers
use integer bisection. Decimal output is for display only. Requires SymPy.
"""
from pathlib import Path
from fractions import Fraction as R
from itertools import combinations
from math import comb
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def require(condition,message):
    if not condition: raise ValueError(message)
z=s.Symbol('z',real=True)
lo,hi=s.Rational(-3,100),s.Rational(-1,250)
q=s.Matrix([1,s.Rational(199,200),s.Rational(2,5),s.Rational(1,400)])
p=s.Rational(759,1000); r0=s.Rational(48959,80000); eta=s.Rational(1,160)
w=s.Matrix([1,-p,-r0-eta*z,z])
a23=(528238451563-2257386450*z)/(200000*(11761-500*z))
a34=2*(1550*z*z-5989657*z+1197903)/(5*(80500*z+48959))
a=[s.Rational(1),s.Rational(28),s.Rational(249,10),a23,s.Rational(4983,200),a34]
lam=-(1550*z+3973343)/200000
pairs=list(combinations(range(4),2))
def generator(rates):
    M=s.zeros(4)
    for (i,j),v in zip(pairs,rates): M[i,j]=v;M[j,i]=v*q[j]/q[i]
    for j in range(4): M[j,j]=-sum(M[:,j])
    return M
M=generator(a); L=M.T
require(all(s.cancel(x)==0 for x in L*w-lam*w),'Eigenvector identity')
require(s.expand((q.T*w)[0])==0,'Bath mean')
require(all(s.cancel(v)==0 for v in s.ones(1,4)*M),'Conservation')
require(M*s.diag(*q)==s.diag(*q)*M.T,'Detailed balance')

# Exact polynomial positivity on [lo, hi]. Transform to t in [0,1].
t=s.Symbol('t')
def bernstein(poly):
    P=s.Poly(s.expand(poly.subs(z,lo+(hi-lo)*t)),t); n=P.degree()
    if n==s.S.NegativeInfinity: return [s.Rational(0)]
    c=[P.nth(j) for j in range(n+1)]
    return [s.factor(sum(c[j]*s.Rational(comb(k,j),comb(n,j)) for j in range(k+1))) for k in range(n+1)]
sign_receipts={}
def positive(expr,label):
    num,den=s.fraction(s.cancel(expr))
    bc,dc=bernstein(num),bernstein(den)
    if all(v<0 for v in dc): bc=[-v for v in bc];dc=[-v for v in dc]
    require(all(v>0 for v in bc+dc),'Continuous sign proof: '+label)
    sign_receipts[label]={'numerator_Bernstein':list(map(str,bc)),'denominator_Bernstein':list(map(str,dc))}
for k,v in enumerate(a): positive(v,'rate_'+str(pairs[k]))
positive(a23-224,'a23_above_224'); positive(a34-10,'a34_above_10')
positive(225-a23,'a23_below_225');positive(12-a34,'a34_below_12')
positive(lam+20,'lambda_above_minus20'); positive(-19-lam,'lambda_below_minus19')
positive(w[2]-w[1],'w3_above_w2');positive(w[3]-w[2],'w4_above_w3');positive(-w[3],'w4_negative')
trip=[]
Fmargin=list(map(s.Rational,['230','1/50','19/10','118']))
Gmargin=list(map(s.Rational,['196','13/100','2/25','3/100']))
for i,j,k in combinations(range(4),3):
    A,B,C=M[i,k],M[j,k],M[i,j];u=q[k]/q[j];tt=q[j]/q[i]
    F=s.factor(A*(1+tt*u/2)-B*(1+u/2)-C*(1-tt));G=s.factor(A*tt-B)
    name=''.join(str(v+1) for v in [i,j,k]);positive(-F,'minusF_'+name);positive(-G,'minusG_'+name)
    idx=len(trip)
    positive(-F-Fmargin[idx],'F_margin_'+name);positive(-G-Gmargin[idx],'G_margin_'+name)
    H=A+B-C*(1+tt);positive(H if idx<2 else -H,'branch_'+name)
    trip.append({'triplet':name,'F':str(F),'G':str(G),'F_upper':str(-Fmargin[idx]),'G_upper':str(-Gmargin[idx])})

# Comparison generator has smaller conductances throughout the whole family.
base_rates=[s.Rational(1),s.Rational(28),s.Rational(249,10),s.Rational(224),s.Rational(4983,200),s.Rational(10)]
Mbase=generator(base_rates)
H=s.diag(*q)*(Mbase.T+50*s.eye(4)); require(H==H.T,'Symmetric inertia matrix')
U=H.copy();pivots=[]
for k in range(4):
    pivot=U[k,k];require(pivot!=0,'Nonzero LDL pivot');pivots.append(pivot)
    for i in range(k+1,4):
        for j in range(k+1,4):U[i,j]=s.factor(U[i,j]-U[i,k]*U[k,j]/pivot)
require(sum(1 for v in pivots if v<0)==2,'Exactly two baseline eigenvalues below -50')
# Dirichlet comparison + minmax then places two family modes below -50;
# the explicit eigenvalue in (-20,-19) must be the simple slow mode.

# Independent algebraic check on the remaining quadratic spectrum.
# Reversibility makes its roots real; Q(-50)>0 and mean<-50 put both below -50.
fast_sum=s.cancel(s.trace(M)-lam)
fast_product=s.cancel((fast_sum**2-(s.trace(M*M)-lam**2))/2)
positive(-fast_sum-100,'fast_roots_mean_below_minus50')
positive(2500+50*fast_sum+fast_product,'fast_quadratic_at_minus50')

# Exact interval primitives.
def rat(v): return R(str(v))
def point(x): return (rat(x),rat(x))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    vs=[x*y for x in a for y in b];return(min(vs),max(vs))
def div(a,b):
    require(b[0]*b[1]>0,'Division interval excludes zero')
    return mul(a,(1/b[1],1/b[0]))
def log_unit(x,N=45):
    require(1<=x<=2,'log range');v=(x-1)/(x+1)
    lower=2*sum((v**(2*k+1)/R(2*k+1) for k in range(N)),R(0))
    remainder=2*v**(2*N+1)/(R(2*N+1)*(1-v*v))
    return(lower,lower+remainder)
def log_bound(x):
    x=rat(x); require(x>=1,'positive logarithm only');k=0
    while x>2:x/=2;k+=1
    return add(log_unit(x),mul(point(k),log_unit(R(2))))
E=[point(0)]+[log_bound(1/rat(qi)) for qi in q[1:]]
e2,e3,e4=E[1:];total=add(add(e2,e3),e4)
# D(0) = A - B*z; D(beta)=Z(beta)^2*b'(beta).
A=add(mul(point(1-p-r0),total),mul(point(4),add(mul(point(p),e2),mul(point(r0),e3))))
B=sub(mul(point(3+eta),e4),add(mul(point(1-eta),e2),mul(point(1+3*eta),e3)))
require(B[0]>0,'Positive threshold denominator')
zc=div(A,B)
require(rat(lo)<zc[0]<zc[1]<rat(hi),'Onset in continuous family')
require(R('-0.017940112')<zc[0] and zc[1]<R('-0.017940111'),'Published onset bracket')
# General single-turn lemma hypotheses: E4>E2+E3 and c13>c23.
require(sub(e4,add(e2,e3))[0]>0,'Separated top energy')
rinterval=add(point(r0),mul(point(eta),(rat(lo),rat(hi))))
c13=mul(add(point(1),rinterval),e3)
c23=mul(sub(point(p),rinterval),sub(e3,e2))
require(c13[0]>c23[1]>0,'Single-turn coefficient dominance')
# Cross-check threshold via direct covariance at z=0 and z=1, avoiding pair formula.
def D0(zz):
    wi=[point(v.subs(z,zz)) for v in w]
    sw=(R(0),R(0));swe=(R(0),R(0))
    for a0,ee in zip(wi,E): sw=add(sw,a0);swe=add(swe,mul(a0,ee))
    return sub(mul(sw,total),mul(point(4),swe))
for zz in [s.Rational(0),s.Rational(1),lo,hi]:
    first=sub(A,mul(B,point(zz)));second=D0(zz)
    require(max(first[0],second[0])<=min(first[1],second[1]),'Two D0 formulas agree')
require(D0(lo)[0]>0 and D0(hi)[1]<0,'Opposite onset sides')

# Derivative at the critical point for the quantitative onset law.
wc=[point(1),point(-p),neg(add(point(r0),mul(point(eta),zc))),zc]
dcrit=point(0)
for i,j in pairs:
    cij=mul(sub(wc[i],wc[j]),sub(E[j],E[i]))
    dcrit=sub(dcrit,mul(cij,add(E[i],E[j])))
require(dcrit[0]>0,'Positive derivative at onset')
width_coefficient=div(B,dcrit)
height_coefficient=div(mul(B,B),mul(point(32),dcrit))

# Predictive threshold for the SAME finite temperature pair across the family.
def nthroot(x,n,digits=24):
    x=rat(x);scale=10**digits;lower,upper=0,scale
    while upper-lower>1:
        mid=(upper+lower)//2
        if mid**n*x.denominator<=x.numerator*scale**n:lower=mid
        else:upper=mid
    require(lower**n*x.denominator<=x.numerator*scale**n<=upper**n*x.denominator,'Integer root enclosure')
    return(R(lower,scale),R(upper,scale))
def affine_overlap(n):
    powers=[nthroot(qi,n) for qi in q]
    Z=(R(0),R(0))
    for v in powers:Z=add(Z,v)
    intercept=sub(powers[0],add(mul(point(p),powers[1]),mul(point(r0),powers[2])))
    slope=sub(powers[3],mul(point(eta),powers[2]))
    return div(intercept,Z),div(slope,Z)
u0,u1=affine_overlap(1000);v0,v1=affine_overlap(100)
gap0=sub(u0,v0);gap1=sub(u1,v1);require(gap1[0]>0,'Finite-pair slope positive')
zpair=div(neg(gap0),gap1)
require(zc[1]<zpair[0]<zpair[1]<rat(hi),'Finite-pair threshold above onset and within family')
# Discrete cases test both sides, including a region with reversal for another pair.
require(add(gap0,mul(gap1,point(hi)))[0]>0,'Reference companion finite witness')
require(add(gap0,mul(gap1,point(lo)))[1]<0,'No fixed-pair witness on lower endpoint')

out={'verified':True,'parameter_interval':[str(lo),str(hi)],'q':list(map(str,q)),'left_mode':list(map(str,w)),'lambda':str(lam),'down_rates':list(map(str,a)),'triplets':trip,'continuous_sign_certificates':sign_receipts,'baseline_LDL_pivots':list(map(str,pivots)),'spectral_conclusion':'lambda2 in (-20,-19); lambda3,lambda4<-50 throughout closed parameter interval','onset_threshold_interval':list(map(str,zc)),'finite_pair_threshold_interval':list(map(str,zpair)),'finite_pair_gap_slope_interval':list(map(str,gap1)),'onset_A_interval':list(map(str,A)),'onset_B_interval':list(map(str,B)),'critical_D_derivative_interval':list(map(str,dcrit)),'peak_beta_linear_coefficient_interval':list(map(str,width_coefficient)),'peak_height_quadratic_coefficient_interval':list(map(str,height_coefficient)),'onset_decimal_display':float(sum(zc)/2),'finite_pair_decimal_display':float(sum(zpair)/2),'theorem_scope':'Ordinary asymptotic cooling iff z>zcrit on specified family; unique cooling amplitude maximum then; no heating reversal on family.'}
(ROOT/'research/mechanism_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('ALL CONTINUOUS-FAMILY EXACT CHECKS PASSED')
print('onset threshold:',out['onset_decimal_display'])
print('fixed finite-pair threshold:',out['finite_pair_decimal_display'])
print('baseline LDL pivots:',list(map(str,pivots)))
print('continuous rational sign certificates:',len(sign_receipts))

print('onset width coefficient:',float(sum(width_coefficient)/2))
print('onset height coefficient:',float(sum(height_coefficient)/2))
