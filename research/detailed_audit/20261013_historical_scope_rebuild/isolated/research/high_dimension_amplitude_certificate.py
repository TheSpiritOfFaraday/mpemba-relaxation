"""Exact clone identities, repair bounds, and an illustrative five-state family member.
Arbitrary-n thermal amplitude limits are analytic in the companion TeX proof.
"""
from pathlib import Path
from math import comb
import sympy as s,json,itertools,time
started=time.time();checks=[];receipts={}
def zero(expr,name):
 assert s.factor(expr)==0,name;checks.append(name)
def require(cond,name):
 assert cond,name;checks.append(name)
e=s.Symbol('e',positive=True);tt=s.Symbol('t');r=s.Symbol('r',positive=True)
p=json.loads(Path(__file__).with_name('quarter_limit_certificate.json').read_text())
M4=s.Matrix([[s.sympify(v,locals={'e':e}) for v in row] for row in p['generator']])
def bernstein_positive(expr,name):
 num,den=s.fraction(s.cancel(expr));data=[]
 for poly in [num,den]:
  poly=s.Poly(s.expand(poly.subs(e,tt/30)),tt);n=poly.degree()
  data.append([s.factor(sum(poly.nth(j)*s.Rational(comb(k,j),comb(n,j)) for j in range(k+1))) for k in range(n+1)])
 if all(v<0 for v in data[1]):data=[[-v for v in row] for row in data]
 require(all(v>0 for row in data for v in row),name)
 receipts[name]={'numerator':list(map(str,data[0])),'denominator':list(map(str,data[1]))}
for i,limit in enumerate([s.Rational(220,169),s.Rational(220,169),s.Rational(132,169)]):
 bernstein_positive(2-M4[i,3],f'base a{i+1}4 <2 over whole interval')
 zero(s.limit(M4[i,3],e,0)-limit,f'bounded rate limit {i+1}')
A,B,C,t,u,rr=s.symbols('A B C t u rr')
F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
zero(F.subs(u,rr*u)-F-(rr-1)*u*G/2,'two-low clone correction identity')
a,q0,qi,rj,rk,mu=s.symbols('a q0 qi rj rk mu',positive=True)
AA=a;BB=mu*rj;CC=a;T=q0*rj/qi;U=rk/rj
zero(AA*T-BB-rj*(a*q0/qi-mu),'one-low clone G identity')
zero(AA*(1+T*U/2)-BB*(1+U/2)-CC*(1-T)-(rj+rk/2)*(a*q0/qi-mu),'one-low clone F identity')
Q=s.Rational(13,5)-e/3+e**2;q=s.Matrix([1,1-e/3,s.Rational(3,5),e**2*r])/Q
kij={(0,1):100,(0,2):10,(1,2):11}
expect=[(-1-98*e/3-e**2*r/2)/Q,(-s.Rational(24,5)-e**2*r)/Q,(-s.Rational(21,5)+3*e-e**2*r/2)/Q]
for (i,j),expr in zip(kij,expect):
 A=q[i]*(i+1);B=q[j]*(j+1);C=q[i]*kij[i,j];t=q[j]/q[i];u=q[3]/q[j]
 F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t)
 zero(F-expr,f'auxiliary low pair {i+1}{j+1} clone F formula')
 bernstein_positive((-F-s.Rational(1,3)).subs(r,0),f'auxiliary low pair {i+1}{j+1} F<-1/3 (worst r=0)')
# Five-state exact member at epsilon=1/1000, clone weights2/3,1/3.
e0=s.Rational(1,1000);base=M4.subs(e,e0);q4=s.Matrix([1,1-e0/3,s.Rational(3,5),e0**2]);q4/=sum(q4)
rclone=[s.Rational(2,3),s.Rational(1,3)];qn=s.Matrix([*q4[:3],q4[3]*rclone[0],q4[3]*rclone[1]])
bar=s.zeros(5);H=s.zeros(5)
for i,j in itertools.combinations(range(5),2):
 if j<3:bar[i,j]=base[i,j];bar[j,i]=base[j,i];H[i,j]=qn[i]*kij[i,j];H[j,i]=qn[j]*kij[i,j]
 elif i<3:bar[i,j]=base[i,3];bar[j,i]=rclone[j-3]*base[3,i];H[i,j]=qn[i]*(i+1);H[j,i]=qn[j]*(i+1)
 else:bar[i,j]=4*rclone[i-3];bar[j,i]=4*rclone[j-3];H[i,j]=4;H[j,i]=4*qn[j]/qn[i]
for mat in [bar,H]:
 for j in range(5):mat[j,j]=-sum(mat[i,j] for i in range(5) if i!=j)
w0=s.Matrix([s.sympify(v,locals={'e':e}).subs(e,e0) for v in p['w']]);w=s.Matrix([*w0[:3],w0[3],w0[3]])
require(bar.T*w==-3*w,'five-state exact lumped left slow mode')
v=s.Matrix([0,0,0,1,-1]);require(bar*v==-(4+sum(base[i,3] for i in range(3)))*v,'five-state exact clone-only fast mode')
M=bar+6*e0**2*H
require(s.ones(1,5)*M==s.zeros(1,5),'five-state conservation')
require(M*s.diag(*qn)==s.diag(*qn)*M.T,'five-state detailed balance')
require(all(M[i,j]>0 for i in range(5) for j in range(5) if i!=j),'five-state all rates positive')
for i,j,k in itertools.combinations(range(5),3):
 A,B,C=M[i,k],M[j,k],M[i,j];t=qn[j]/qn[i];u=qn[k]/qn[j]
 require(A*t-B<0,f'five-state strict G{i+1}{j+1}{k+1}')
 require(A*(1+t*u/2)-B*(1+u/2)-C*(1-t)<0,f'five-state strict F{i+1}{j+1}{k+1}')
Cn=214;gap=s.Rational(1,3)-12*Cn*e0**2;require(gap>s.Rational(1,4),'five-state perturbed gap lower bound >1/4')
out={'verified':True,'checks':checks,'bernstein_certificates':receipts,'exact_check_count':len(checks),'elapsed_seconds':time.time()-started,'five_state_epsilon':str(e0),'five_state_gap_lower_bound':str(gap),'scope':'Generic cloning identities; rigorous rate/repair bounds over 0<epsilon<=1/30; exact five-state strict positive family member. Thermal advantage limits proved analytically in high_dimension_amplitude.tex.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['checks','bernstein_certificates']},indent=2))
