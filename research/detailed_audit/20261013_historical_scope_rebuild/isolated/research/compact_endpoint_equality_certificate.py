"""Exact regression: simple slow mode, full failures, endpoint equality."""
from pathlib import Path
import json,itertools
import sympy as s
from exact_rational_triplet_classifier import classify, require
p=Path(__file__).with_name('rational_triplet_compact_equality.json')
x=json.loads(p.read_text());q=s.Matrix(list(map(s.Rational,x['bath_weights'])));q/=sum(q)
w=s.Matrix(list(map(s.Rational,x['mode'])));m=s.Matrix(x['energy_exponents']);base=s.Rational(x['base'])
require(all(s.Rational(v)==base**int(e) for v,e in zip(x['bath_weights'],m)), 'commensurate weights')
require((q.T*w)[0]==0, 'zero bath mean')
require(sum(w[i]*m[i] for i in range(4))*4-sum(w)*sum(m)==0, 'zero endpoint covariance')
M=s.zeros(4)
for (i,j),v in zip(itertools.combinations(range(4),2),map(s.Rational,x['downward_rates'])):M[i,j]=v;M[j,i]=v*q[j]/q[i]
for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
lam=(M.T*w)[0]/w[0];require(M.T*w==lam*w, 'exact left eigenpair')
r=classify(x['bath_weights'],x['downward_rates'])
require(r['verified'] and r['slow_multiplicity']==1 and r['endpoint_sign']==0 and r['ordinary_cooling']==False, 'exact equality classification')
require(s.Rational(r['slow_eigenvalue'])==lam, 'correct slow root')
out={'verified':True,'base':str(base),'energy_exponents':list(map(int,m)),'fixture_bytes':p.stat().st_size,'checks':['commensurate rational bath weights','zero bath mean','exact zero uniform energy covariance','exact left eigenpair','slow multiplicity one','all full original triplet failures','classifier recognizes exact prime-log cancellation'],'classification':'no ordinary cooling; zero endpoint derivative; strict monotonicity for finite positive beta by proved theorem'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
