"""Symbolic inverse construction for an asymptotic amplitude family (exploration)."""
import sympy as s,json,itertools
from pathlib import Path
e=s.Symbol('e',positive=True);q=s.Matrix([1,1-e/3,s.Rational(3,5),e*e]);q/=sum(q);D=s.diag(*q)
w=s.Matrix([1,-s.Rational(4,5),-s.Rational(1,3)-4*e/9+5*e**3/3,-e])
U=s.Matrix.hstack(*(q.T.col_join((D*w).T)).nullspace())
x,y=s.symbols('x y');B=s.Matrix([[1+x,y],[y,1-x]])/2;M=-D*U*B*U.T
rows=[]
for i,j,k in itertools.combinations(range(4),3):
 A,B0,C=M[i,k],M[j,k],M[i,j];t=q[j]/q[i];u=q[k]/q[j]
 rows.extend([s.factor(A*(1+t*u/2)-B0*(1+u/2)-C*(1-t)),s.factor(A*t-B0)])
sol=s.solve([rows[2],rows[7]],[x,y]);B=s.simplify(B.subs(sol));M=s.simplify(M.subs(sol))
print('U',U);print('coordinates',{str(k):str(s.factor(v)) for k,v in sol.items()});print('B',B)
print('detB',s.factor(B.det()))
for n,row in enumerate(rows):print('row',n,s.factor(row.subs(sol)))
out=dict(U=[[str(s.factor(v)) for v in U.row(i)] for i in range(4)],B=[[str(s.factor(v)) for v in B.row(i)] for i in range(2)],rows=[str(s.factor(v.subs(sol))) for v in rows],detB=str(s.factor(B.det())))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
