"""Exact boundary regressions, including a missed thermal derivative root."""
from pathlib import Path
import json
import sympy as s
from exact_polynomial_reference import sign_changes_open
from exact_commensurate_classifier import OrderedField, X

checks=[]
def require(ok,label):
    if not ok: raise ValueError(label)
    checks.append(label)

x=X;half=s.Rational(1,2);one=s.Rational(1)
Z=1+x+x*x+x**3
N=x+(-s.Rational(8,17)+s.Rational(1,10**30))*x**3
H=s.Poly(s.expand(s.diff(N,x)*Z-N*s.diff(Z,x)),x,domain=s.QQ)
require(H.eval(half)>0 and H.eval(one)<0,'Opposite exact endpoint signs imply an interior sign change')
intervals=H.intervals(eps=s.Rational(1,10**12))
inside=[(a,b) for (a,b),m in intervals if half<a and b<one]
require(not inside and H.eval((half+one)/2)<0,
        'Old finite-width reference would omit the near-bath root and report only negative values')
require(sign_changes_open(H,half,one)==1,'New reference counts the near-bath sign change exactly')
for name,expr,left,right,expected in [
 ('excluded endpoints',(x-half)*(x-one),half,one,0),
 ('even stationary tangency',(x-s.Rational(3,4))**2,half,one,0),
 ('odd repeated stationary root',(x-s.Rational(3,4))**3,half,one,1),
 ('two simple interior roots',(x-s.Rational(2,3))*(x-s.Rational(3,4)),half,one,2),
 ('constant nonzero',s.Integer(-7),half,one,0),
 ('near zero endpoint',x-s.Rational(1,10**40),s.Rational(0),half,1),
 ('near infinite temperature endpoint',x-1+s.Rational(1,10**40),half,one,1),
]:
 require(sign_changes_open(s.Poly(expr,x,domain=s.QQ),left,right)==expected,name)
field=OrderedField(s.Rational(1))
require(field.root_count_open(field.odd_part(H),half,one)==1,
        'Main exact classifier independently handles the newly exposed boundary case')
alpha=s.CRootOf(s.Symbol('z')**3-s.Symbol('z')-1,0)
af=OrderedField(alpha)
P=s.Poly((x-half)*(x-alpha/2)**2*(x-s.Rational(3,4))**3,x,domain=af.K)
require(af.root_count_open(af.odd_part(P),half,one)==1,
        'Algebraic coefficient field retains odd multiplicity and excludes endpoint')
require(af.root_count_open(P,half,one)==2,
        'Algebraic strong-root count includes even interior multiplicity once')
out={'passed':True,'checks':checks,'count':len(checks),
 'finding':'Supporting classifier_audit.py formerly used finite-width isolating boxes to decide whether roots were interior. A root just above x=1/2 could be skipped. Production exact_commensurate_classifier.py was unaffected.',
 'repair':'Replace that reference with exact rational odd-factor root counts and exact endpoint-factor removal.'}
Path('research/detailed_audit/polynomial_reference_repair.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'checks':len(checks)},indent=2))
