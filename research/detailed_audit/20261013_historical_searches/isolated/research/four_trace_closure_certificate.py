"""Exact counterexample: tracing the highest state need not preserve full F/G."""
import sympy as s,itertools,json
from pathlib import Path
q=s.Matrix([151,150,100,99])/500;delta=s.Rational(1,40);pairs=list(itertools.combinations(range(4),2))
k=[1+delta,s.Rational(401,201)+delta,delta,2+delta,delta,2+delta]
# A tiny all-downward-one repair makes every original inequality strict.
repair=s.Rational(1,10**6)
K=s.zeros(4)
for (i,j),a in zip(pairs,k):K[i,j]=K[j,i]=a+repair/q[i]
def tests(K,q):
 out=[]
 for i,j,k in itertools.combinations(range(len(q)),3):
  F=(q[i]+q[k]/2)*K[i,k]-(q[j]+q[k]/2)*K[j,k]-(q[i]-q[j])*K[i,j]
  G=q[j]*(K[i,k]-K[j,k]);out.append(dict(triplet=[i+1,j+1,k+1],F=F,G=G))
 return out
original=tests(K,q);assert all(z['F']<0 and z['G']<0 for z in original)
d=sum(q[i]*K[i,3] for i in range(3));Kt=K[:3,:3]+q[3]*K[:3,3]*K[3,:3]/d
traced=tests(Kt,q[:3]);assert traced[0]['F']>0
# Diagonal of reduced K has no meaning; only offdiagonal rates are used.
out=dict(verified=True,bath_weights=list(map(str,q)),reduced_conductances=[str(K[i,j]) for i,j in pairs],downward_rates=[str(q[i]*K[i,j]) for i,j in pairs],original_tests=[{k:str(v) if k in ['F','G'] else v for k,v in t.items()} for t in original],trace_tests=[{k:str(v) if k in ['F','G'] else v for k,v in t.items()} for t in traced],trace_F_float=float(traced[0]['F']),construction='k=(1+d,401/201+d,d,2+d,d,2+d), d=1/40, plus constant-downward1 generator times1/10^6',scope='Elimination of highest state before subtracting minimal reset does not preserve all F, even for strict original failure. No assertion about the reset-subtracted trace.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
