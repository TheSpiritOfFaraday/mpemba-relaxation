"""Exact four-state example of a purely transient chi-squared reversal.
All rates, mode projections, and polynomial signs are rational. Temperatures
are beta_h=1/3, beta_c=2/3, beta_b=1. Crossing times are logarithms of two
algebraic roots. This example is unrestricted, not in the full triplet cone.
"""
import sympy as s,json,itertools
from pathlib import Path
from exact_commensurate_classifier import require
x=s.symbols('x');q=s.Matrix([s.Rational(1,8)**i for i in range(4)]);q/=sum(q);D=s.diag(*q)
ph=s.Matrix([s.Rational(1,2)**i for i in range(4)]);ph/=sum(ph)
pc=s.Matrix([s.Rational(1,4)**i for i in range(4)]);pc/=sum(pc)
N=s.expand((x-s.Rational(1,2))*(x-s.Rational(1,4))*(x-s.Rational(1,8)))
base=s.Matrix([N.coeff(x,i) for i in range(4)]);energy=s.Matrix(range(4));energy-=s.ones(4,1)*(q.T*energy)[0]
w=base+energy/1000
require((q.T*w)[0]==0,'bath orthogonality')
v=s.Matrix.hstack(*(q.T.col_join((D*w).T).col_join(ph.T)).nullspace())
require(v.shape==(4,1),'middle direction dimension')
u=s.Matrix.hstack(*(q.T.col_join((D*w).T).col_join((D*v).T)).nullspace())
require(u.shape==(4,1),'fast direction dimension')
vec=[w,v,u];proj=[z*z.T*D/(z.T*D*z)[0] for z in vec]
require(sum(proj,s.zeros(4))==s.eye(4)-s.ones(4,1)*q.T,'complete weighted orthogonal projections')
core=(proj[1]+2*proj[2]).T
rho=1+max([s.Rational(0)]+[core[i,j]/q[i] for i in range(4) for j in range(4) if i!=j])
M=rho*(q*s.ones(1,4)-s.eye(4))-core
require(s.ones(1,4)*M==s.zeros(1,4) and M*D==D*M.T,'generator identities')
require(all(M[i,j]>0 for i in range(4) for j in range(4) if i!=j),'strictly positive rates')
for j,z in enumerate(vec):require(M.T*z==(-rho-j)*z,'ordered exact eigenpair')
abc=[s.factor(((z.T*ph)[0]**2-(z.T*pc)[0]**2)/(z.T*D*z)[0]) for z in vec]
A,B,C=abc;F=A+B*x+C*x*x;xm=-B/(2*C);disc=B*B-4*A*C
require(A>0 and B<0 and C>0,'normal eventual ordering and alternating coefficients')
require(A+B+C>0,'normal initial ordering')
require(0<xm<1 and F.subs(x,xm)<0,'strict transient reversal')
require(disc>0,'two simple crossings')
root1=(-B-s.sqrt(disc))/(2*C);root2=(-B+s.sqrt(disc))/(2*C)
# Rational Sturm counts certify both roots lie between zero and one, without
# evaluating the radicals numerically for acceptance.
require(s.Poly(F,x).count_roots(0,1)==2,'two roots in physical time interval')
trip=[]
for i,j,k in itertools.combinations(range(4),3):
 aa,bb,cc=M[i,k],M[j,k],M[i,j];t=q[j]/q[i];uu=q[k]/q[j]
 FF=aa*(1+t*uu/2)-bb*(1+uu/2)-cc*(1-t);GG=aa*t-bb;kap=(aa+bb)/(cc*(1+t))
 trip.append(dict(states=[i+1,j+1,k+1],F=str(FF),G=str(GG),predicate=bool((kap>=1 and FF>0) or (kap<=1 and GG>0))))
out=dict(verified=True,q=list(map(str,q)),energies='(0,3 log2,6 log2,9 log2)',beta_hot='1/3',beta_cold='2/3',beta_bath='1',
 generator=[[str(v) for v in M.row(i)] for i in range(4)],downward_rates=[str(M[i,j]) for i,j in itertools.combinations(range(4),2)],
 eigenvalues=['0',str(-rho),str(-rho-1),str(-rho-2)],chi_squared_coefficients=list(map(str,abc)),
 scaled_difference='A+B*x+C*x^2, x=exp(-2t)',critical_x=str(xm),critical_polynomial_value=str(s.factor(F.subs(x,xm))),
 crossing_times=['-log((-B+sqrt(B^2-4AC))/(2C))/2','-log((-B-sqrt(B^2-4AC))/(2C))/2'],
 displayed_crossing_times=[str(s.N(-s.log(root2)/2,15)),str(s.N(-s.log(root1)/2,15))],
 displayed_coefficients=[str(s.N(v,12)) for v in abc],triplets=trip,
 scope='Exact purely transient weighted-Euclidean reversal with normal eventual ordering; not asserted to satisfy triplet failure')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['generator','downward_rates','triplets','chi_squared_coefficients','critical_polynomial_value']},indent=2))
