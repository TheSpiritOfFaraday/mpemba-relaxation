"""Independent regression audit; does not modify either classifier."""
import itertools,json,random,sys,time
from pathlib import Path
import sympy as s
import mpmath as mp
from exact_rational_triplet_classifier import classify,log_prime_bounds,algebraic_bounds
from exact_commensurate_classifier import classify as classify_old,OrderedField
start=time.time(); cases=[]; PAIRS=list(itertools.combinations(range(4),2))
def run(name,q,r,expected=None,cross=None):
 out=classify(q,r)
 assert out['verified'],name
 if expected is not None: assert out['ordinary_cooling']==expected,name
 if cross:
  old=classify_old(cross[0],cross[1],r)
  for k in ['ordinary_cooling','ordinary_heating','slow_multiplicity']:assert out[k]==old[k],(name,k)
  assert not old['strong_cooling_count'] and not old['strong_heating_count']
 cases.append({'name':name,'ordinary_cooling':out['ordinary_cooling'],'endpoint_sign':out.get('endpoint_sign'),'slow_multiplicity':out['slow_multiplicity'],'crosschecked':bool(cross),'interval_rounds':out.get('interval_rounds')})
 return out
q=list(map(s.Rational,['1','199/200','2/5','1/400']));r=list(map(s.Rational,['1','28','249/10','1799/8','4983/200','251/25']))
a=run('primary rational certified cooling',q,r,True)
for scale in [s.Rational(7,13),s.Rational(1000),s.Rational(1,104729)]:
 b=run('common weight scale '+str(scale),[z*scale for z in q],r,True)
 assert b['nonzero_prime_log_coefficients']==a['nonzero_prime_log_coefficients']
for scale in [s.Rational(7,13),s.Rational(1000),s.Rational(1,104729)]:run('time scale '+str(scale),q,[z*scale for z in r],True)
for base,m in [(s.Rational(1,2),[0,1,2,3]),(s.Rational(2,3),[0,1,3,7])]:
 q0=[base**k for k in m];q=[z/sum(q0) for z in q0]
 run('reset '+str(base),q,[q[i] for i,j in PAIRS],False,(base,m))
 for k in range(3):
  r=[q[i]+(s.Rational(7,3)*q[i]/(q[i]+q[j]) if (i,j)==(k,k+1) else 0) for i,j in PAIRS]
  out=run('accelerated adjacent pair '+str((base,k)),q,r,False,(base,m))
  assert out['slow_multiplicity']==2
 run('constant downward '+str(base),q,[1]*6,False,(base,m))
 w=s.Matrix([1,1,1,-sum(q[:3])/q[3]]);D=s.diag(*q);norm=(w.T*D*w)[0]
 M=s.Matrix(q)*s.ones(1,4)-s.eye(4)+s.Rational(1,10)*D*w*w.T/norm
 run('simple slow with repeated faster roots '+str(base),q,[M[i,j] for i,j in PAIRS],False,(base,m))
 # Small arbitrary perturbation of strict full-premise constant-downward generator.
 rng=random.Random(81279)
 for k in range(4):
  r=[1+s.Rational(rng.randrange(-5,6),1000) for _ in range(6)]
  run('strict perturbed full-premise '+str((base,k)),q,r,False,(base,m))
# Algebraic enclosure with nontrivial embeddings and nonmonotone coefficient signs.
x=s.Symbol('z')
for idx in range(3):
 alpha=s.CRootOf(x**3-3*x+1,idx);f=OrderedField(7*alpha/13)
 for v in [alpha-1,2-alpha**2,3*alpha**2-5*alpha+1]:
  exact=f.convert(v)
  for _ in range(12):f.interval=f.interval.refine()
  lo,hi=algebraic_bounds(f,exact); numerical=s.N(v,80)
  assert s.N(lo,80)<=numerical<=s.N(hi,80)
  assert f.sign(exact)==s.sign(numerical)
cases.append({'name':'three real embeddings and scaled generator intervals','passed':True})
mp.mp.dps=90
for p in [2,3,5,199,104729,2147483647]:
 last=None
 for terms in [4,8,16,32]:
  lo,hi=log_prime_bounds(p,terms)
  val=mp.log(p); conv=lambda z:mp.mpf(int(z.p))/int(z.q)
  assert conv(lo)<val<conv(hi)
  if last:assert last[0]<=lo<=hi<=last[1]
  last=(lo,hi)
cases.append({'name':'six prime logs with nested exact enclosures and independent 90 digit check','passed':True})
invalid=[('float weights',[1.,.5,.25,.125],[1]*6),('equal weights',[1,1,s.Rational(1,4),s.Rational(1,8)],[1]*6),('zero rate',[1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)],[0,1,1,1,1,1]),('outside premise',[1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)],[1,100,1,1,1,1])]
for name,q,r in invalid:
 try:classify(q,r)
 except ValueError:cases.append({'name':'reject '+name,'passed':True})
 else:raise AssertionError(name)
result={'verified':True,'case_count':len(cases),'elapsed_seconds':time.time()-start,'cases':cases,'limitation':'Equality and bounded-round inconclusive path are not exercised here; parent is constructing exact equality.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
