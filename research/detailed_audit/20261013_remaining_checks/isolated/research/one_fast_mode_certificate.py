"""Exact symbolic audit of the arbitrary-n one-fast-mode classification.
Finite dimension checks are supporting evidence; the support induction
in one_fast_mode_classification.tex proves the arbitrary-n assertion.
"""
from pathlib import Path
import itertools,json,time
import sympy as s
started=time.time();checks=[]
def zero(v,name):
    assert s.factor(v)==0,name;checks.append(name)
def matrixzero(A,name):
    assert all(s.factor(v)==0 for v in A),name;checks.append(name)
def tests(M,q,i,j,k):
    A,B,C=M[i,k],M[j,k],M[i,j]
    return A*(1+q[k]/(2*q[i]))-B*(1+q[k]/(2*q[j]))-C*(1-q[j]/q[i]),A*q[j]/q[i]-B
p,q,r=s.symbols('p q r',positive=True);u,v,w=s.symbols('u v w',real=True)
rho,delta,N=s.symbols('rho delta N',positive=True)
A=rho*p-delta*p*u*w/N;B=rho*q-delta*q*v*w/N;C=rho*p-delta*p*u*v/N
zero(A*q/p-B+delta*q*w*(u-v)/N,'generic rank-one G')
zero(A*(1+r/(2*p))-B*(1+r/(2*q))-C*(1-q/p)-delta*((p-q)*u*v-w*(p*u-q*v+r*(u-v)/2))/N,'generic rank-one F')
for n in range(4,9):
    weights=s.symbols(f'q1:{n+1}',positive=True);D=s.diag(*weights);Q=sum(weights)
    R=s.Matrix(weights)*s.ones(1,n)-Q*s.eye(n)
    Bdiff=s.zeros(n-1,n);H=s.zeros(n,n-1)
    for i in range(n-1):Bdiff[i,i]=-1;Bdiff[i,i+1]=1
    for i in range(n):
        for j in range(i):H[i,j]=1
    matrixzero(R*R+Q*R,f'n={n} reset projector identity')
    for k in range(n-1):
        z=s.zeros(n,1);z[k]=-1/weights[k];z[k+1]=1/weights[k+1]
        norm=(z.T*D*z)[0];J=-D*z*z.T/norm;M=rho*R+delta*J
        matrixzero(J*J+J,f'n={n} pair={k+1} rank-one eigenvalue identity')
        zero(s.trace(J)+1,f'n={n} pair={k+1} rank-one trace')
        matrixzero(R*J+Q*J,f'n={n} pair={k+1} reset action left')
        matrixzero(J*R+Q*J,f'n={n} pair={k+1} reset action right')
        matrixzero(s.ones(1,n)*M,f'n={n} pair={k+1} conservation')
        matrixzero(M*D-D*M.T,f'n={n} pair={k+1} detailed balance')
        expected=-rho*Q*s.eye(n-1);expected[k,k]-=delta
        if k>0:expected[k-1,k]+=delta*weights[k+1]/(weights[k]+weights[k+1])
        if k<n-2:expected[k+1,k]+=delta*weights[k]/(weights[k]+weights[k+1])
        matrixzero(Bdiff*M.T*H-expected,f'n={n} pair={k+1} difference generator')
        d=delta/(weights[k]+weights[k+1])
        for i,j,l in itertools.combinations(range(n),3):
            F,G=tests(M,weights,i,j,l);Fe=Ge=0
            if (i,j)==(k,k+1):Fe=-(weights[k]-weights[k+1])*d
            if (j,l)==(k,k+1):Fe=-(weights[k]+weights[k+1]/2)*d;Ge=-weights[k]*d
            zero(F-Fe,f'n={n} pair={k+1} F{i+1},{j+1},{l+1}')
            zero(G-Ge,f'n={n} pair={k+1} G{i+1},{j+1},{l+1}')
# Sharp dimension caveat: all three entries of the fast mode nonzero.
q=s.Matrix([s.Rational(3,5),s.Rational(3,10),s.Rational(1,10)]);z=s.Matrix([1,2,-12]);D=s.diag(*q)
M=q*s.ones(1,3)-s.eye(3)-D*z*z.T/(z.T*D*z)[0]
assert all(M[i,j]>0 for i in range(3) for j in range(3) if i!=j)
checks.append('n=3 caveat all rates positive')
F,G=tests(M,q,0,1,2);zero(F,'n=3 caveat F=0');zero(G+s.Rational(2,9),'n=3 caveat G=-2/9')
x=s.Symbol('x');zero(M.charpoly(x).as_expr()-x*(x+1)*(x+2),'n=3 caveat spectrum')
result={'verified':True,'exact_check_count':len(checks),'elapsed_seconds':time.time()-started,'checks':checks,'scope':'Generic rank-one identities; symbolic pair matrices n=4..8; n=3 necessity caveat. Arbitrary-n proof in companion TeX.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
