"""Find fixed-mode PSD dual obstruction for a rational n-state test mode."""
import itertools,json
from pathlib import Path
import sympy as s
import numpy as np
from scipy.optimize import linprog
q=s.Matrix([1,s.Rational(1,2),s.Rational(1,4),s.Rational(1,8),s.Rational(1,16)]);n=len(q);q/=sum(q);D=s.diag(*q);Q=s.Rational(0);H=s.zeros(n,n-1)
for j in range(n-1):
 Q+=q[j]
 for i in range(j+1):H[i,j]=1
 H[j+1,j]=-Q/q[j+1]
v=s.Matrix([1,s.Rational(1,28),1,1]);w=H*v # threshold b=1/14, so violates by factor two
assert w[0]+w[1]-2*w[2]<0
U=s.Matrix.hstack(*(q.T.col_join((D*w).T)).nullspace());d=n-2
pairs=list(itertools.combinations(range(n),2));triples=list(itertools.combinations(range(n),3));bpairs=list(itertools.combinations_with_replacement(range(d),2));Bs=[]
for i,j in bpairs:
 B=s.zeros(d);B[i,j]=B[j,i]=1;Bs.append(B)
rows=[];names=[]
for i,j,k in triples:
 row=[];grow=[]
 for B in Bs:
  M=-D*U*B*U.T;A,B0,C=M[i,k],M[j,k],M[i,j]
  row.append(A*(1+q[k]/(2*q[i]))-B0*(1+q[k]/(2*q[j]))-C*(1-q[j]/q[i]));grow.append(A*q[j]/q[i]-B0)
 rows.extend([row,grow]);names.extend(['F'+''.join(str(z+1) for z in (i,j,k)),'G'+''.join(str(z+1) for z in (i,j,k))])
C=s.Matrix(rows);target=s.Matrix([int(i==j) for i,j in bpairs]);sol=linprog(np.ones(len(rows)),A_eq=np.array(C.T,dtype=float),b_eq=np.array(target,dtype=float).ravel(),bounds=[(0,None)]*len(rows),method='highs')
out=dict(feasible_dual=sol.success,q=list(map(str,q)),w=list(map(str,w)),arithmetic_margin=str(w[0]+w[1]-2*w[2]))
if sol.success:
 support=np.flatnonzero(sol.x>1e-7).tolist();eq=C.T[:,support];values=list(s.linsolve((eq,target)));assert len(values)==1
 a=s.Matrix(values[0]);assert all(x>=0 for x in a);assert eq*a==target
 out.update(verified=True,dual={names[i]:str(a[j]) for j,i in enumerate(support)},meaning='Positive combination of original F/G constraints equals trace(B); full premise incompatible with B positive definite for this mode.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
