"""Exact boundary, multiplicity, algebraic-sign and model tests."""
import json,time
from pathlib import Path
import sympy as s
from exact_commensurate_classifier import classify,OrderedField,PAIRS,X,require
start=time.time();tests=[]
def make_mode(w,base=s.Rational(1,2),m=(0,1,2,3)):
 q=s.Matrix([base**i for i in m]);q/=sum(q);D=s.diag(*q);w=s.Matrix(w)
 w-=s.ones(4,1)*(q.T*w)[0];norm=(w.T*D*w)[0]
 eps=1/(10*(1+max(abs(w[i]*w[j])/norm for i,j in PAIRS)))
 M=q*s.ones(1,4)-s.eye(4)+eps*D*w*w.T/norm
 require(all(M[i,j]>0 for i,j in PAIRS),'positive test rates')
 return [M[i,j] for i,j in PAIRS]
def check(name,base,m,rates,expected):
 result=classify(base,m,rates)
 for k,v in expected.items():require(result[k]==v,name+': '+k)
 tests.append(dict(name=name,passed=True,result=result));return result
no=dict(ordinary_cooling=False,ordinary_heating=False,strong_cooling_count=0,strong_heating_count=0)
check('both strong directions','1/2',[0,1,2,3],make_mode([-s.Rational(3,32),s.Rational(11,16),-s.Rational(3,2),1]),dict(ordinary_cooling=True,ordinary_heating=True,strong_cooling_count=1,strong_heating_count=1,slow_multiplicity=1))
N=s.expand((X-s.Rational(1,2))*(X-s.Rational(3,4))**2)
check('tangent strong root','1/2',[0,1,2,3],make_mode([N.coeff(X,i) for i in range(4)]),dict(ordinary_cooling=True,ordinary_heating=False,strong_cooling_count=1,strong_heating_count=0,cooling_sign_changing_stationary_points=2))
N=s.expand((X-s.Rational(1,2))*(X-1)*(X-s.Rational(1,4)))
check('infinite-temperature root excluded','1/2',[0,1,2,3],make_mode([N.coeff(X,i) for i in range(4)]),dict(strong_cooling_count=0,strong_heating_count=1))
w2,w3,w4=s.symbols('w2 w3 w4');N=1+w2*X+w3*X**2+w4*X**3;Z=1+X+X**2+X**3;der=s.diff(N,X)*Z-N*s.diff(Z,X)
sol=s.solve([N.subs(X,s.Rational(1,2)),der.subs(X,s.Rational(3,4)),s.diff(der,X).subs(X,s.Rational(3,4))],[w2,w3,w4])
check('stationary tangency is not an effect','1/2',[0,1,2,3],make_mode([1,sol[w2],sol[w3],sol[w4]]),no)
check('bath threshold equality','1/3',[0,1,2,3],make_mode([1,-3,0,0],s.Rational(1,3)),no)
q=s.Matrix([1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)]);q/=sum(q);D=s.diag(*q)
prep=s.Matrix([s.Rational(3,4)**i for i in range(4)]);prep/=sum(prep)
U=s.Matrix.hstack(s.ones(4,1),D.inv()*prep);P=s.eye(4)-U*(U.T*D*U).inv()*U.T*D
L=-s.eye(4)+s.ones(4,1)*q.T+P/1000;M=L.T
check('double slow mode','1/2',[0,1,2,3],[M[i,j] for i,j in PAIRS],dict(slow_multiplicity=2,strong_cooling_count=1,strong_heating_count=0,ordinary_cooling=True))
check('triple slow reset','1/2',[0,1,2,3],[q[i] for i,j in PAIRS],dict(no,slow_multiplicity=3))
check('irrational cubic slow mode','1/2',[0,1,2,3],[1,2,3,4,5,6],dict(no,slow_multiplicity=1))
check('zero rates connected chain','1/2',[0,1,2,3],[1,0,0,1,0,1],dict(no,slow_multiplicity=1))
check('nonuniform commensurate spacing','1/2',[0,1,3,7],[1,2,3,4,5,6],dict(slow_multiplicity=1))
check('scaled algebraic slow root','1/2',[0,3,5,8],['49/5','3/4','89/100','61/10','37/5',82],dict(slow_multiplicity=1))
# AlgebraicField.is_positive is not a real-embedding oracle: alpha-2 is
# negative despite a positive leading coefficient in its representation.
t=s.symbols('t');alpha=s.CRootOf(t**3-t-1,0);field=OrderedField(alpha)
require(field.sign(field.convert(alpha-2))==-1,'ordered embedding regression')
tests.append(dict(name='ordered algebraic embedding regression',passed=True))
field=OrderedField(3*alpha)
require(field.sign(field.convert(3*alpha-4))==-1,'scaled real embedding regression')
require(field.sign(field.convert(3*alpha-3))==1,'scaled real embedding positive sign')
tests.append(dict(name='scaled algebraic embedding regression',passed=True))
for name,base,m,rates in [('noninteger exponent','1/2',[0,1,2.5,3],[1]*6),('disconnected','1/2',[0,1,2,3],[1,0,0,0,0,1]),('equal energies','1/2',[0,1,1,3],[1]*6)]:
 try:classify(base,m,rates)
 except ValueError:tests.append(dict(name=name,passed=True))
 else:raise ValueError(name+' should be rejected')
out=dict(verified=True,test_count=len(tests),tests=tests,elapsed_seconds=time.time()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(verified=True,test_count=len(tests),names=[t['name'] for t in tests],elapsed_seconds=time.time()-start),indent=2))
