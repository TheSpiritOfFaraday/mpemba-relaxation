import sympy as s,numpy as np
from scipy.optimize import linprog
from itertools import combinations
from exact_mode_realizability import PAIRS,TRIPLES
q=s.Matrix(list(map(s.Rational,['5','4','2','1'])));w=s.Matrix([162,-18,71,-880]);D=s.diag(*q);U=s.Matrix.hstack(*(q.T.col_join((D*w).T)).nullspace());print('U',U)
bases=[s.Matrix([[1,0],[0,0]]),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,0],[0,1]])]
down=s.Matrix([[-q[i]*(U[i,:]*B*U[j,:].T)[0] for B in bases] for i,j in PAIRS]);rows=[]
for i,j,k in TRIPLES:
 t=q[j]/q[i];u=q[k]/q[j];A=down[PAIRS.index((i,k)),:];B=down[PAIRS.index((j,k)),:];C=down[PAIRS.index((i,j)),:]
 rows.extend([A*(1+t*u/2)-B*(1+u/2)-C*(1-t),A*t-B])
R=s.Matrix.vstack(*rows)
for ids in combinations(range(8),3):
 V=R[list(ids),:].T
 if V.det()==0:continue
 z=V.inv()*s.Matrix([-1,0,-1])
 if min(z)>=0:print('DUAL',ids,z.T)
print('R',R)
