"""Exact symbolic identities and rational realization for the classification.

This checks algebra and a sharpness example. The all-energy zero-count and
classification proofs are in universal_classification.tex, not inferred here
from finite samples. Uses SymPy rational arithmetic; no discovery output read.
"""
import json
from pathlib import Path
import sympy as s

def require(condition, label):
    if not condition:
        raise ValueError(label)

a,b,c=s.symbols('a b c',positive=True)
f2,f3,f4=s.symbols('f2 f3 f4',positive=True)
f=s.Matrix([1,f2,f3,f4]);Z=sum(f);energies=[0,a,b,c]
def diff_beta(expr):
    return s.expand(sum(-energies[i]*f[i]*s.diff(expr,f[i]) for i in range(1,4)))
pi=f/Z
g=[s.cancel(diff_beta(pi[i])) for i in range(4)]
u=s.cancel(g[1]/g[0]);v=s.cancel(g[2]/g[0])
slope=s.cancel(diff_beta(v)/diff_beta(u))
expected=(b*c*(c-b)*f4/f2-a*b*(b-a))/(a*b*(b-a)+a*c*(c-a)*f4/f3)
require(s.cancel(slope-expected)==0,'slope identity')
require(s.cancel(u-(-a+(b-a)*f3+(c-a)*f4)/(a+b*f3/f2+c*f4/f2))==0,'u identity')
require(s.cancel(v-(f3/f2)*(-b+(a-b)*f2+(c-b)*f4)/(a+b*f3/f2+c*f4/f2))==0,'v identity')

q=s.Matrix([1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8)])
p=q/sum(q);D=s.diag(*p)
w=s.Matrix([-s.Rational(3,32),s.Rational(11,16),-s.Rational(3,2),1])
N=(w.T*D*w)[0];eps=s.Rational(1,100)
M=-s.eye(4)+p*s.ones(1,4)+eps*D*w*w.T/N
require((p.T*w)[0]==0,'bath orthogonality')
require(s.ones(1,4)*M==s.zeros(1,4),'column sums')
require(M*D==D*M.T,'detailed balance')
require(all(M[i,j]>0 for i in range(4) for j in range(4) if i!=j),'positive rates')
require(M.T*w==(-1+eps)*w,'slow eigenvector')
lam=s.Symbol('lambda')
require(s.factor(M.charpoly(lam).as_expr())==lam*(lam+1)**2*(100*lam+99)/100,'full spectrum')
x=s.Symbol('x');poly=sum(w[i]*x**i for i in range(4))
require(s.expand(poly-(x-s.Rational(1,4))*(x-s.Rational(1,2))*(x-s.Rational(3,4)))==0,'three exact Gibbs zeros')
H=s.Poly(s.expand(s.diff(poly,x)*(1+x+x*x+x**3)-poly*s.diff(1+x+x*x+x**3,x)),x)
require(H.count_roots(0,s.Rational(1,2))==1,'one heating stationary point')
require(H.count_roots(s.Rational(1,2),1)==1,'one cooling stationary point')
out={'verified':True,'scope':'symbolic coordinate and slope identities; rational sharpness realization',
     'q':list(map(str,q)),'w':list(map(str,w)),'weighted_norm_squared':str(N),
     'epsilon':str(eps),'generator':[[str(M[i,j]) for j in range(4)] for i in range(4)],
     'spectrum':['0','-99/100','-1','-1'],'Gibbs_numerator':str(s.factor(poly)),
     'derivative_polynomial':str(H.as_expr()),'stationary_roots_cooling':1,'stationary_roots_heating':1,
     'strong_beta':['log(4/3)/log(2)','2'],'triplet_constraints_claimed':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
