"""Exact continuous-parameter certificate for a 25-percent-limit family.
All inequalities cover 0<e<=1/30 by rational Bernstein coefficients.
The limiting thermal argument is recorded separately in the report fragment.
"""
import sympy as s,json,itertools,time
from math import comb
from pathlib import Path
start=time.time();e=s.Symbol('e',positive=True);tt=s.Symbol('t')
q=s.Matrix([1,1-e/3,s.Rational(3,5),e*e]);q/=sum(q);D=s.diag(*q)
w=s.Matrix([1,-s.Rational(4,5),-s.Rational(1,3)-4*e/9+5*e**3/3,-e])
def require(ok,msg):
 if not ok:raise ValueError(msg)
U=s.Matrix.hstack(*(q.T.col_join((D*w).T)).nullspace())
x,y=s.symbols('x y');B=s.Matrix([[1+x,y],[y,1-x]])/2;M0=-D*U*B*U.T
triples=list(itertools.combinations(range(4),3));rows=[]
for i,j,k in triples:
 A,B0,C=M0[i,k],M0[j,k],M0[i,j];t=q[j]/q[i];u=q[k]/q[j]
 rows.extend([s.factor(A*(1+t*u/2)-B0*(1+u/2)-C*(1-t)),s.factor(A*t-B0)])
linear_system=s.linear_eq_to_matrix([rows[2],rows[7]],[x,y])
sol=s.solve([rows[2],rows[7]],[x,y]);B=s.simplify(B.subs(sol));M=s.simplify(-D*U*(B/e**2)*U.T+3*(q*s.ones(1,4)-s.eye(4)))
receipts={}
def coeff(poly):
 P=s.Poly(s.expand(poly.subs(e,tt/30)),tt);n=P.degree()
 return [s.factor(sum(P.nth(j)*s.Rational(comb(k,j),comb(n,j)) for j in range(k+1))) for k in range(n+1)]
def positive(expr,name):
 num,den=s.fraction(s.cancel(expr));P=s.Poly(num,e);Q=s.Poly(den,e)
 en=min(m[0] for m,c in P.terms());ed=min(m[0] for m,c in Q.terms())
 num=s.cancel(num/e**en);den=s.cancel(den/e**ed);nn,dd=coeff(num),coeff(den)
 if all(v<0 for v in dd):nn=[-v for v in nn];dd=[-v for v in dd]
 require(all(v>0 for v in nn+dd),'positive Bernstein coefficients: '+name)
 receipts[name]=dict(e_power=en-ed,numerator_coefficients=list(map(str,nn)),denominator_coefficients=list(map(str,dd)))
require(s.cancel((q.T*w)[0])==0,'bath orthogonality')
require(all(s.cancel(v)==0 for v in M.T*w+3*w),'exact slow mode')
require(all(s.cancel(v)==0 for v in s.ones(1,4)*M),'column conservation')
require(all(s.cancel(v)==0 for v in M*D-D*M.T),'detailed balance')
positive(B[0,0],'B11');positive(B.det(),'detB')
linear_det=s.factor(linear_system[0].det())
positive(linear_det if linear_det.subs(e,s.Rational(1,30))>0 else -linear_det,'linear_system_nonsingular')
for i in range(4):
 for j in range(4):
  if i!=j:positive(M[i,j],'rate_%d%d'%(i+1,j+1))
trip=[]
for i,j,k in triples:
 A,B0,C=M[i,k],M[j,k],M[i,j];t=q[j]/q[i];u=q[k]/q[j]
 F=s.factor(A*(1+t*u/2)-B0*(1+u/2)-C*(1-t));G=s.factor(A*t-B0)
 for name,v in [('F',F),('G',G)]:
  if v!=0:positive(-v,'minus%s%d%d%d'%(name,i+1,j+1,k+1))
 trip.append(dict(triplet=[i+1,j+1,k+1],F=str(F),G=str(G)))
gram=U.T*D*U;fast_trace=s.cancel(s.trace(B*gram)/e**2);fast_det=s.cancel(B.det()*gram.det()/e**4);gap=s.cancel(fast_det/fast_trace)
positive(gap-s.Rational(1,3),'fast_gap_above_one_third')
# A fixed strict interior generator permits an arbitrarily small perturbation.
kk=[100,10,1,11,2,3];H=s.zeros(4)
for (i,j),k in zip(itertools.combinations(range(4),2),kk):H[i,j]=q[i]*k;H[j,i]=q[j]*k
for j in range(4):H[j,j]=-sum(H[i,j] for i in range(4) if i!=j)
for i,j,k in triples:
 A,B0,C=H[i,k],H[j,k],H[i,j];t=q[j]/q[i];u=q[k]/q[j]
 F=s.factor(A*(1+t*u/2)-B0*(1+u/2)-C*(1-t));G=s.factor(A*t-B0)
 positive(-F,'interior_minusF%d%d%d'%(i+1,j+1,k+1));positive(-G,'interior_minusG%d%d%d'%(i+1,j+1,k+1))
limitM=M.applyfunc(lambda v:s.limit(v,e,0));out=dict(verified=True,interval='0<e<=1/30',q_raw=['1','1-e/3','3/5','e^2'],w=list(map(str,w)),
 U=[[str(s.factor(v)) for v in U.row(i)] for i in range(4)],B=[[str(s.factor(v)) for v in B.row(i)] for i in range(2)],
 slow_eigenvalue='-3',fast_gap_lower_bound='1/3',fast_trace=str(fast_trace),fast_det=str(fast_det),
 generator=[[str(s.factor(v)) for v in M.row(i)] for i in range(4)],triplets=trip,
 limit_generator=[[str(v) for v in limitM.row(i)] for i in range(4)],limit_spectrum={str(k):v for k,v in limitM.eigenvals().items()},
 certificates=receipts,certificate_count=len(receipts),elapsed_seconds=time.time()-start,
 scope='Exact continuous positivity/triplet/spectral proof; thermal25percent limit requires accompanying analytic argument')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['generator','certificates','B','U','triplets']},indent=2))
