"""Exact identities for the full-triplet four-state slow-degeneracy theorem.

Analytic sign cases are proved in full_triplet_degenerate.tex. Unnormalized
positive q are allowed below; then reset eigenvalue is -rho*sum(q).
"""
from pathlib import Path
import itertools,json,time
import sympy as s
started=time.time(); checks=[]
def zero(expr,name):
    assert s.factor(expr)==0,name
    checks.append(name)
def matrixzero(expr,name):
    assert all(s.factor(v)==0 for v in expr),name
    checks.append(name)
q=s.symbols('q1:5',positive=True); z=s.symbols('z1:5',real=True)
rho,delta=s.symbols('rho delta',positive=True); lam=s.Symbol('lambda')
D=s.diag(*q); Q=sum(q);R=s.Matrix(q)*s.ones(1,4)-Q*s.eye(4)
N=sum(qi*zi**2 for qi,zi in zip(q,z))
M=rho*R-delta*D*s.Matrix(z)*s.Matrix(z).T/N
triplets=list(itertools.combinations(range(4),3))
def tests(M,i,j,k):
    A,B,C=M[i,k],M[j,k],M[i,j]
    return s.factor(A*(1+q[k]/(2*q[i]))-B*(1+q[k]/(2*q[j]))-C*(1-q[j]/q[i])),s.factor(A*q[j]/q[i]-B)
for i,j,k in triplets:
    F,G=tests(M,i,j,k)
    zero(G+delta*q[j]*z[k]*(z[i]-z[j])/N,f'rank-one G identity {i+1}{j+1}{k+1}')
    Fexpected=delta*((q[i]-q[j])*z[i]*z[j]-z[k]*(q[i]*z[i]-q[j]*z[j]+q[k]*(z[i]-z[j])/2))/N
    zero(F-Fexpected,f'rank-one F identity {i+1}{j+1}{k+1}')
# Difference operator and right inverse f1=0.
B=s.Matrix([[-1,1,0,0],[0,-1,1,0],[0,0,-1,1]])
H=s.Matrix([[0,0,0],[1,0,0],[1,1,0],[1,1,1]])
for k in range(3):
    zz=s.zeros(4,1);zz[k]=-1/q[k];zz[k+1]=1/q[k+1]
    nn=(zz.T*D*zz)[0]; J=-D*zz*zz.T/nn; MM=rho*R+delta*J
    zero(sum(q[i]*zz[i] for i in range(4)),f'pair {k+1} zero bath mean')
    matrixzero(s.ones(1,4)*MM,f'pair {k+1} column conservation')
    matrixzero(MM*D-D*MM.T,f'pair {k+1} detailed balance')
    zero(J[k,k+1]-q[k]/(q[k]+q[k+1]),f'pair {k+1} downward rate')
    zero(J[k+1,k]-q[k+1]/(q[k]+q[k+1]),f'pair {k+1} upward rate')
    zero(MM.charpoly(lam).as_expr()-lam*(lam+rho*Q)**2*(lam+rho*Q+delta),f'pair {k+1} exact spectrum')
    expectedK=-rho*Q*s.eye(3);expectedK[k,k]-=delta
    if k>0:expectedK[k-1,k]+=delta*q[k+1]/(q[k]+q[k+1])
    if k<2:expectedK[k+1,k]+=delta*q[k]/(q[k]+q[k+1])
    matrixzero(B*MM.T*H-expectedK,f'pair {k+1} Metzler difference matrix')
    for i,j,l in triplets:
        F,G=tests(MM,i,j,l);Fe=Ge=s.Integer(0)
        if k==0 and (i,j)==(0,1):Fe=-delta*(q[0]-q[1])/(q[0]+q[1])
        if k==1:
            if (i,j,l)==(0,1,2):Fe=-delta*(2*q[1]+q[2])/(2*(q[1]+q[2]));Ge=-delta*q[1]/(q[1]+q[2])
            if (i,j,l)==(1,2,3):Fe=-delta*(q[1]-q[2])/(q[1]+q[2])
        if k==2 and (j,l)==(2,3):Fe=-delta*(2*q[2]+q[3])/(2*(q[2]+q[3]));Ge=-delta*q[2]/(q[2]+q[3])
        zero(F-Fe,f'pair {k+1} converse F{i+1}{j+1}{l+1}')
        zero(G-Ge,f'pair {k+1} converse G{i+1}{j+1}{l+1}')
result={'verified':True,'exact_checks':checks,'elapsed_seconds':time.time()-started,'scope':'Symbolic rank-one tests, converse, conservation, reversibility, spectrum and stochastic monotonicity; analytic support-case proof in companion TeX.','assumptions':'q1>q2>q3>q4>0; rho>0; delta>0. Normalized theorem uses sum(q)=1.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
