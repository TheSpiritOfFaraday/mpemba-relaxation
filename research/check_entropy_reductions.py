"""Exact local identities supporting the analytic entropy/dimensional review.
This checks algebra, not the universal analytic theorems or the open conjecture.
"""
from pathlib import Path
import itertools,json
import sympy as s
checks=[]
def check(name,value):
    if not value: raise AssertionError(name)
    checks.append(name)
def zero(expr): return s.factor(expr)==0
# A 3-coordinate generator's exterior-square operator: derive by differentiating minors.
b=s.Matrix(3,3,s.symbols('b:9')); u=s.Matrix(s.symbols('u:3'));v=s.Matrix(s.symbols('v:3'))
w=s.Matrix([u[i]*v[j]-u[j]*v[i] for i,j in [(0,1),(0,2),(1,2)]])
C=s.Matrix([[b[0,0]+b[1,1],b[1,2],-b[0,2]], [b[2,1],b[0,0]+b[2,2],b[0,1]],[-b[2,0],b[1,0],b[1,1]+b[2,2]]])
dw=s.Matrix([(b*u)[i]*v[j]+u[i]*(b*v)[j]-(b*u)[j]*v[i]-u[j]*(b*v)[i] for i,j in [(0,1),(0,2),(1,2)]])
check('Exterior-square differential operator has exactly the stated crossing signs',all(zero(z) for z in dw-C*w))
Ct=C.subs({b[0,2]:0,b[2,0]:0})
check('Tridiagonal specialization removes both negative off-diagonal crossing terms',all(Ct[i,j] in [0,b[0,1],b[1,0],b[1,2],b[2,1]] for i in range(3) for j in range(3) if i!=j))
a,beta,mu=s.symbols('a b mu',positive=True)
check('Angle lower bound follows from an exact square',zero(mu**2*(a+beta)**2-4*a*beta*((a+beta)*mu-a*beta)-((a+beta)*mu-2*a*beta)**2))
z=s.symbols('z',real=True)
check('Elementary KL threshold gives 1-z^2-z^2/4=1-5z^2/4',zero(1-z*z-z*z/4-(1-5*z*z/4)))
# General symbolic compression with two rare-tail states. F uses original physical rates.
p,q,r,h,j=s.symbols('p q r h j',positive=True);weights=[p,q,r,h,j]
K=s.zeros(5)
for ii,jj in itertools.combinations(range(5),2): K[ii,jj]=K[jj,ii]=s.Symbol('k%d%d'%(ii,jj))
M=s.diag(*weights)*K
for ii in range(5):M[ii,ii]=-sum(M[jj,ii] for jj in range(5) if jj!=ii)
Q=h+j;W=[p,q,r,Q];T=s.zeros(4)
for ii in range(3):
 for jj in range(3):
  if ii!=jj:T[ii,jj]=M[ii,jj]
 T[ii,3]=(h*M[ii,3]+j*M[ii,4])/Q
 T[3,ii]=M[3,ii]+M[4,ii]
for ii in range(4):T[ii,ii]=-sum(T[jj,ii] for jj in range(4) if jj!=ii)
def fg(N,w,i,j,k):
 A,B,C=N[i,k],N[j,k],N[i,j]
 return (A*(1+w[k]/(2*w[i]))-B*(1+w[k]/(2*w[j]))-C*(1-w[j]/w[i]),A*w[j]/w[i]-B)
for ii,jj in itertools.combinations(range(3),2):
 F,G=fg(T,W,ii,jj,3);orig=[fg(M,weights,ii,jj,k) for k in [3,4]]
 check('Compressed G weighted-average identity %d%d'%(ii,jj),zero(G-(h*orig[0][1]+j*orig[1][1])/Q))
 correction=sum(weights[k]*(Q-weights[k])*orig[k-3][1]/(2*Q*weights[jj]) for k in [3,4])
 check('Compressed F weighted-average plus nonpositive correction %d%d'%(ii,jj),zero(F-(h*orig[0][0]+j*orig[1][0])/Q-correction))
check('Compressed conservation and detailed balance',all(zero(v) for v in T*s.diag(*W)-s.diag(*W)*T.T) and all(zero(v) for v in s.ones(1,4)*T))
# Symbolic weighted graph identity with three contrast coordinates.
x=s.symbols('x:3',positive=True);eta=s.symbols('eta:3',positive=True);y=s.Matrix(s.symbols('y:3'));ell=s.symbols('ell:3',nonnegative=True)
B=s.zeros(3);energy=0
for (i,j),e in zip(itertools.combinations(range(3),2),ell):
 B[i,j]=e/eta[i];B[j,i]=e/eta[j]
 B[i,i]-=e*x[j]/(eta[i]*x[i]);B[j,j]-=e*x[i]/(eta[j]*x[j])
 energy+=e*x[i]*x[j]*(y[i]/x[i]-y[j]/x[j])**2
check('Graph quadratic form is an exact sum of edge squares',zero(-(y.T*s.diag(*eta)*B*y)[0]-energy))
check('Graph kernel contains prescribed mode exactly',all(zero(t) for t in B*s.Matrix(x)))
# First-three arithmetic in unnormalized Helmert coordinates, tails cancel.
xx,yy,tail=s.symbols('xx yy tail');ww=[xx+yy+tail,-p*xx/q+yy+tail,-(p+q)*yy/r+tail]
check('Arithmetic violation is the displayed strict coordinate inequality',zero(ww[0]+ww[1]-2*ww[2]-(2*(p+q+r)/r)*(yy-(p-q)*r*xx/(2*q*(p+q+r)))))
out={'passed':True,'count':len(checks),'checks':checks,'scope':__doc__.strip()}
Path('research/detailed_audit/entropy_reduction_identities.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
