"""Exact boundary and sharpness checks for the uniform cooling threshold."""
import sympy as s,json
from pathlib import Path
def require(ok,msg):
    if not ok:raise ValueError(msg)
x=s.Symbol('x');Z=1+x+x*x+x**3
N_boundary=1-3*x;H_boundary=s.factor(s.diff(N_boundary,x)*Z-N_boundary*s.diff(Z,x))
require(s.expand(H_boundary-2*(x-1)*(3*x*x+3*x+2))==0,'threshold tangency identity')
q=s.Matrix([1,s.Rational(1,4),s.Rational(1,16),s.Rational(1,64)]);p=q/sum(q);D=s.diag(*p)
w=s.Matrix([1,-4,0,0]);eta=(w.T*D*w)[0];eps=s.Rational(1,100)
M=-s.eye(4)+p*s.ones(1,4)+eps*D*w*w.T/eta
require(all(M[i,j]>0 for i in range(4) for j in range(4) if i!=j),'positive rates')
require(s.ones(1,4)*M==s.zeros(1,4) and M*D==D*M.T,'Markov reversibility')
require(M.T*w==-s.Rational(99,100)*w,'slow eigenvector')
l=s.Symbol('lambda');require(s.expand(M.charpoly(l).as_expr()-l*(l+1)**2*(l+s.Rational(99,100)))==0,'simple slow spectrum')
for i,j,k in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:require(M[i,k]/p[i]==M[j,k]/p[j],'G equality')
N=1-4*x;xh=s.Rational(31,32);xc=s.Rational(15,16)
bh=(N/Z).subs(x,xh);bc=(N/Z).subs(x,xc)
require(bh<0 and bc<bh,'strict cooling witness')
out=dict(verified=True,energies='(0,log(4),2log(4),3log(4))',bath_beta='1',
    q=list(map(str,q)),left_mode=list(map(str,w)),spectrum=['0','-99/100','-1','-1'],
    generator=[[str(M[i,j]) for j in range(4)] for i in range(4)],
    hot_beta='log(32/31)/log(4)',cold_beta='log(16/15)/log(4)',
    hot_overlap=str(bh),cold_overlap=str(bc),strict_amplitude_margin=str(bh-bc),
    boundary_derivative_polynomial=str(H_boundary),G_conditions='all equal zero',
    original_F_conditions_claimed=False,scope='Sharpness in the G-only family and exact equality-boundary polynomial')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
