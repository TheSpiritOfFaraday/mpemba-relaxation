"""Exact rational projection of an unaccepted floating inverse-search candidate."""
import json,itertools
from pathlib import Path
import sympy as s
import mpmath as mp
mp.mp.dps=80
p=json.load(open('research/general_n_inverse_candidate.json'));q=s.Matrix([s.Rational(str(x)) for x in p['q']]);q/=sum(q);n=len(q);w=s.Matrix([s.Rational(str(x)) for x in p['w']]);w-=s.ones(n,1)*(q.T*w)[0];D=s.diag(*q);norm=(w.T*D*w)[0];P=s.eye(n)-s.ones(n,1)*q.T-w*(w.T*D)/norm
conv=lambda z:mp.mpf(int(z.p))/int(z.q)
U=s.Matrix([[s.Rational(mp.nstr(mp.mpf(str(p['U'][i][j]))/mp.sqrt(conv(q[i])),40)) for j in range(n-2)] for i in range(n)]);U=P*U;B=s.Matrix([[s.Rational(str(v)) for v in row] for row in p['B']]);assert all(B[i,i]>sum(abs(B[i,j]) for j in range(n-2) if j!=i) for i in range(n-2));M=-D*U*B*U.T
assert M.T*w==s.zeros(n,1) and s.ones(1,n)*M==s.zeros(1,n)
evals=[]
for i,j,k in itertools.combinations(range(n),3):
 F=M[i,k]*(1+q[k]/(2*q[i]))-M[j,k]*(1+q[k]/(2*q[j]))-M[i,j]*(1-q[j]/q[i]);G=M[i,k]*q[j]/q[i]-M[j,k]
 for typ,z in [('F',F),('G',G)]:
  if z>0:evals.append(dict(test=typ+''.join(str(t+1) for t in (i,j,k)),value_approx=str(s.N(z,30))))
assert evals
out=dict(accepted=False,reason='Exact rational projected realization violates original triplet constraints; floating feasibility residual cannot be accepted.',positive_tests=evals,arithmetic_margin=str(w[0]+w[1]-2*w[2]),positive_definite_fast_matrix=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
