"""Exact ordinary/strong asymptotic classification of four-state systems.

Input: distinct integer energy exponents m, rational 0<base<1, six rational
downward rates in 12,13,14,23,24,34 order. E_i=-m_i log(base), beta_b=1.
Rates may vanish if the graph remains connected. All slow multiplicities are
handled. No floating arithmetic is used for acceptance or classification.

For algebraic eigenvalues, coefficient signs use rational interval refinement
at the selected real embedding, NOT AlgebraicField.is_positive (which does not
in general implement that ordered embedding). SymPy 1.14.0 supplies exact
polynomial arithmetic and CRootOf rational isolating intervals.
"""
import argparse,itertools,json
from pathlib import Path
import sympy as s

X,T=s.symbols('x lambda')
PAIRS=list(itertools.combinations(range(4),2))
def require(ok,msg):
 if not ok:raise ValueError(msg)

class OrderedField:
 def __init__(self,alpha):
  self.alpha=alpha;self.rational=bool(alpha.is_Rational)
  # all_roots(radicals=False) may return a rational multiple of CRootOf,
  # rather than a bare CRootOf. Use the root atom as the field generator,
  # so the coefficient polynomial and isolating interval refer to the same
  # algebraic number. This also handles rational affine root expressions.
  if self.rational:
   self.K=s.QQ;self.interval=None
  else:
   atoms=list(alpha.atoms(s.CRootOf))
   require(len(atoms)==1,'one real algebraic root embedding')
   self.generator=atoms[0]
   self.K=s.QQ.algebraic_field(self.generator)
   self.K.from_sympy(alpha)
   self.interval=self.generator._get_interval()
  self.refinements=0
 def convert(self,v):return self.K.from_sympy(s.sympify(v))
 def sign(self,v):
  if not v:return 0
  if self.rational:return 1 if v>0 else -1
  coeff=[s.Rational(c.numerator,c.denominator) for c in v.to_list()]
  while True:
   a=s.Rational(self.interval.a.numerator,self.interval.a.denominator)
   b=s.Rational(self.interval.b.numerator,self.interval.b.denominator)
   lo=hi=s.Rational(0)
   for c in coeff:
    products=[lo*a,lo*b,hi*a,hi*b];lo=min(products)+c;hi=max(products)+c
   if lo>0:return 1
   if hi<0:return -1
   self.interval=self.interval.refine();self.refinements+=1
 def evaluate(self,P,x):return P.rep.eval(self.convert(x))
 def root_count_open(self,P,a,b):
  require(a<b,'ordered root-count interval')
  if P.degree()<=0:return 0
  P=P.sqf_part()
  # Remove endpoint factors explicitly; neither infinite nor bath endpoints
  # are allowed to masquerade as finite nonbath preparations.
  for end in [a,b]:
   while P.degree()>0 and not self.evaluate(P,end):
    P,rem=P.div(s.Poly(X-end,X,domain=self.K));require(rem.is_zero,'endpoint division')
  if P.degree()<=0:return 0
  seq=P.sturm()
  def changes(end):
   signs=[self.sign(self.evaluate(f,end)) for f in seq]
   signs=[z for z in signs if z]
   return sum(a!=b for a,b in zip(signs,signs[1:]))
  result=changes(a)-changes(b);require(result>=0,'Sturm root count')
  return result
 def odd_part(self,P):
  require(not P.is_zero,'nonconstant thermal amplitude')
  result=s.Poly(1,X,domain=self.K)
  for factor,mult in P.sqf_list()[1]:
   if mult%2:result*=factor
  return result

def classify(base,exponents,rates):
 base=s.Rational(base);raw_m=list(map(s.sympify,exponents));rates=list(map(s.Rational,rates))
 require(all(v.is_Integer for v in raw_m),'integer energy exponents without truncation')
 m=list(map(int,raw_m))
 require(0<base<1,'rational bath base lies between zero and one')
 require(len(m)==4 and all(m[i]<m[i+1] for i in range(3)),'four ordered integer exponents')
 m=[v-m[0] for v in m]
 require(len(rates)==6 and all(v>=0 for v in rates),'six nonnegative rates')
 edges=[pair for pair,a in zip(PAIRS,rates) if a>0];seen={0}
 for _ in range(4):
  for i,j in edges:
   if i in seen or j in seen:seen|={i,j}
 require(len(seen)==4,'irreducible transition graph')
 q=s.Matrix([base**v for v in m]);q/=sum(q);D=s.diag(*q);M=s.zeros(4)
 for (i,j),a in zip(PAIRS,rates):M[i,j]=a;M[j,i]=a*q[j]/q[i]
 for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
 require(M*D==D*M.T and s.ones(1,4)*M==s.zeros(1,4),'reversibility and conservation')
 char=s.Poly(M.charpoly(T).as_expr()/T,T);roots=char.all_roots(radicals=False)
 require(all(r.is_real for r in roots),'real reversible spectrum')
 alpha=roots[-1];r=sum(z==alpha for z in roots);field=OrderedField(alpha);K=field.K
 require(field.sign(field.convert(alpha))<0,'negative slow eigenvalue')
 Z=s.Poly(sum(X**v for v in m),X,domain=K);L=M.T
 if r==1:
  # Adjugate columns give a slow left eigenvector. Form them over rational
  # polynomials first and reduce at the selected real algebraic eigenvalue.
  adj=(L-T*s.eye(4)).adjugate();w=None
  for j in range(4):
   trial=[field.convert(adj[i,j].subs(T,alpha)) for i in range(4)]
   if any(trial):w=trial;break
  require(w is not None,'nonzero simple-mode adjugate')
  require(all(sum(field.convert(L[i,j])*w[j] for j in range(4))==field.convert(alpha)*w[i] for i in range(4)),'exact algebraic eigenpair')
  require(not sum(field.convert(q[i])*w[i] for i in range(4)),'zero bath overlap')
  N=s.Poly.from_dict({(mi,):wi for mi,wi in zip(m,w)},X,domain=K)
  derivative=N.diff()*Z-N*Z.diff();thermal_power=1
 else:
  # A repeated root of this rational cubic is rational; the remaining root
  # is rational too. Polynomial spectral projections avoid square roots of q.
  require(alpha.is_Rational,'repeated slow root is rational')
  if r==2:
   fast=roots[0];P=L*(L-fast*s.eye(4))/(alpha*(alpha-fast))
  else:P=s.eye(4)-s.ones(4,1)*q.T
  require(P*P==P and P*s.ones(4,1)==s.zeros(4,1),'exact slow projection')
  H=P*D.inv();require(H==H.T,'symmetric amplitude form')
  f=s.Matrix([X**v for v in m]);N=s.Poly((f.T*H*f)[0],X,domain=K)
  derivative=N.diff()*Z-2*N*Z.diff();thermal_power=2
 odd=field.odd_part(derivative)
 cooling_turns=field.root_count_open(odd,base,s.Rational(1))
 heating_turns=field.root_count_open(odd,s.Rational(0),base)
 strong_cooling=field.root_count_open(N,base,s.Rational(1))
 strong_heating=field.root_count_open(N,s.Rational(0),base)
 return dict(verified=True,method='ordered algebraic field and exact Sturm root counts',
  base=str(base),energy_exponents=m,downward_rates=list(map(str,rates)),
  slow_eigenvalue=str(alpha),slow_multiplicity=r,nonzero_characteristic_polynomial=str(char.as_expr()),
  ordinary_cooling=bool(cooling_turns),ordinary_heating=bool(heating_turns),
  strong_cooling_count=strong_cooling,strong_heating_count=strong_heating,
  cooling_sign_changing_stationary_points=cooling_turns,heating_sign_changing_stationary_points=heating_turns,
  thermal_numerator=str(N.as_expr()),partition_power=thermal_power,
  derivative_numerator=str(derivative.as_expr()),algebraic_sign_refinements=field.refinements,
  scope='Exact for rational rates and commensurate integer energy gaps; asymptotic KL and chi-squared for all multiplicities, any fixed norm for simple modes. No transient classification.')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('input',help='JSON containing base, energy_exponents, downward_rates')
 parser.add_argument('--output');args=parser.parse_args();data=json.loads(Path(args.input).read_text())
 out=classify(data['base'],data['energy_exponents'],data['downward_rates']);render=json.dumps(out,indent=2)+'\n'
 if args.output:Path(args.output).write_text(render)
 print(render,end='')
