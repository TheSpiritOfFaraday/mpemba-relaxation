"""Exact five-state obstruction to a proposed minimum-reset/trace proof route.
This is NOT a counterexample to the first-three slow-mode inequality.
"""
import sympy as s,itertools,json
from pathlib import Path
q=s.Matrix([100,4,3,2,1])/110;pairs=list(itertools.combinations(range(5),2));K=s.zeros(5)
values=[1,s.Rational(999,1000),s.Rational(999,1000),0,1,s.Rational(999,1000),0,2,0,2]
eta=s.Rational(1,10**6)
for (i,j),v in zip(pairs,values):K[i,j]=K[j,i]=v+eta/q[i]
minimum=min(K[i,j] for i,j in pairs)
assert minimum==eta/q[0]
for i,j in pairs:K[i,j]-=minimum;K[j,i]-=minimum
assert min(K[i,j] for i,j in pairs)==0
# Physical input adds reset1; its minimum reduced rate is exactly1.
Kphysical=K.copy()
for i,j in pairs:Kphysical[i,j]+=1;Kphysical[j,i]+=1
assert min(Kphysical[i,j] for i,j in pairs)==1
assert all(Kphysical[i,j]>0 for i,j in pairs)
def tests(K,q):
 out=[]
 for i,j,k in itertools.combinations(range(len(q)),3):
  F=s.factor((q[i]+q[k]/2)*K[i,k]-(q[j]+q[k]/2)*K[j,k]-(q[i]-q[j])*K[i,j]);G=s.factor(q[j]*(K[i,k]-K[j,k]))
  out.append({'triplet':[i+1,j+1,k+1],'F':F,'G':G})
 return out
original=tests(Kphysical,q);normalized=tests(K,q)
assert all(t['F']<0 and t['G']<0 for t in original)
assert original==normalized
Kt=K.copy()
for n in [4,3]:
 d=sum(q[i]*Kt[i,n] for i in range(n));assert d>0
 Kt=Kt[:n,:n]+q[n]*Kt[:n,n]*Kt[n,:n]/d
trace=tests(Kt,q[:3]);assert trace[0]['F']>0
out={'verified':True,'q':list(map(str,q)),'base_reduced_rates':list(map(str,values)),'strict_downward_repair':str(eta),'minimum_subtracted_from_repaired_base':str(minimum),'physical_input_reduced_rates':[str(Kphysical[i,j]) for i,j in pairs],'physical_minimum_reset':'1','all_original_tests_strict':True,'trace_F123':str(trace[0]['F']),'trace_F123_float':float(trace[0]['F']),'scope':'Even after subtracting the minimum reset from a strictly positive full-F/G generator, the trace onto first3 may violate F123. This rules out the proposed zero-energy Schur sufficient lemma; it does not refute the all-n slow-mode arithmetic conjecture.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
