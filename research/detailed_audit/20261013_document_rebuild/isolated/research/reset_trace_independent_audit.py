"""Independent direct-generator/block-Schur verification of trace obstruction."""
from pathlib import Path
import json,itertools
import sympy as s
r=json.loads(Path('research/reset_trace_obstruction.json').read_text())
q=s.Matrix([s.Rational(z) for z in r['q']]);n=5
K=s.zeros(n)
for (i,j),v in zip(itertools.combinations(range(n),2),r['physical_input_reduced_rates']):K[i,j]=K[j,i]=s.Rational(v)
def generator(K):
 M=s.zeros(n)
 for i,j in itertools.combinations(range(n),2):M[i,j]=q[i]*K[i,j];M[j,i]=q[j]*K[i,j]
 for j in range(n):M[j,j]=-sum(M[i,j] for i in range(n) if i!=j)
 return M
M=generator(K);D=s.diag(*q)
assert M*D==D*M.T and s.ones(1,n)*M==s.zeros(1,n)
assert min(M[i,j] for i in range(n) for j in range(n) if i!=j)>0
assert min(K[i,j] for i,j in itertools.combinations(range(n),2))==1
def fg(M,q):
 out=[]
 for i,j,k in itertools.combinations(range(len(q)),3):
  F=M[i,k]*(1+q[k]/(2*q[i]))-M[j,k]*(1+q[k]/(2*q[j]))-M[i,j]*(1-q[j]/q[i])
  G=M[i,k]*q[j]/q[i]-M[j,k]
  out.append((s.factor(F),s.factor(G)))
 return out
original=fg(M,q);assert all(F<0 and G<0 for F,G in original)
R=q*s.ones(1,n)-s.eye(n);N=M-R
assert fg(N,q)==original
assert min(N[i,j] for i in range(n) for j in range(n) if i!=j)==0
assert all(N[i,j]>=0 for i in range(n) for j in range(n) if i!=j)
T=N[:3,:3]-N[:3,3:]*N[3:,3:].inv()*N[3:,:3]
assert s.ones(1,3)*T==s.zeros(1,3) and T*s.diag(*q[:3])==s.diag(*q[:3])*T.T
F,G=fg(T,q[:3])[0];assert F==s.Rational(r['trace_F123'])>0
H=s.zeros(n,n-1)
for k in range(n-1):
 for i in range(k+1):H[i,k]=1
 H[k+1,k]=-sum(q[:k+1])/q[k+1]
eta=H.T*D*H;B=eta.inv()*H.T*D*N.T*H
Bt=B[:2,:2]-B[:2,2:]*B[2:,2:].inv()*B[2:,:2]
H3=H[:3,:2];eta3=H3.T*s.diag(*q[:3])*H3
Beff=eta3.inv()*H3.T*s.diag(*q[:3])*T.T*H3
assert Bt==Beff
z=(q[0]-q[1])*q[2]/(2*q[1]*sum(q[:3]))
C=s.factor(Bt[1,0]+z*(Bt[1,1]-Bt[0,0])-z*z*Bt[0,1]);assert C<0
out={'verified':True,'original_max_F':str(max(f for f,g in original)),'original_max_G':str(max(g for f,g in original)),'tail_block_determinant':str(N[3:,3:].det()),'trace_F':str(F),'trace_G':str(G),'zero_energy_schur_arithmetic_quadratic':str(C),'trace_equals_helmert_zero_schur':True}
Path('research/reset_trace_independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
