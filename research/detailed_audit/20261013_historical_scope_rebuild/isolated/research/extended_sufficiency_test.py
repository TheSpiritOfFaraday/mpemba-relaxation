import sympy as s,json
from exact_mode_realizability import exact_realize
q=s.Matrix([1,s.Rational(4,5),s.Rational(2,5),s.Rational(1,5)])
H=s.Matrix([[1,1,1],[-q[0]/q[1],1,1],[0,-(q[0]+q[1])/q[2],1],[0,0,-sum(q[:3])/q[3]]])
for b in [s.Rational(23,1000),s.Rational(1,40),s.Rational(1,20),s.Rational(1,10),s.Rational(1,2),1,2,10]:
 for c in [s.Rational(1,1000000),s.Rational(1,1000),s.Rational(1,100),s.Rational(1,10),1,10,100]:
  w=H*s.Matrix([1,b,c]);res=exact_realize(q,w)
  if not res['feasible']:
   out=dict(q=list(map(str,q)),b=str(b),c=str(c),w=list(map(str,w)),arith=str(w[0]+w[1]-2*w[2]),weighted_mean=str((q.T*w)[0]),result=res)
   print(json.dumps(out,indent=2));json.dump(out,open('research/extended_sufficiency_counterexample.json','w'),indent=2);raise SystemExit
print('all feasible')
