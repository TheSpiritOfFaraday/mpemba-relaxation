"""Recompute stored borderline LP search mode and repair rational roundoff."""
import json,itertools
from pathlib import Path
import sympy as s
import mpmath as mp
mp.mp.dps=100
p=json.load(open('research/general_n_arithmetic_search_result.json'))['best'];q=s.Matrix([s.Rational(str(x)) for x in p['q']]);q/=sum(q);n=len(q);pairs=list(itertools.combinations(range(n),2));k=list(map(lambda z:s.Rational(str(z)),p['reduced_rates']));rho=1-min(k);M=s.zeros(n);H=s.zeros(n)
for (i,j),v in zip(pairs,k):
 M[i,j]=q[i]*(v+rho);M[j,i]=q[j]*(v+rho)
 H[i,j]=1;H[j,i]=q[j]/q[i]
for j in range(n):M[j,j]=-sum(M[i,j] for i in range(n) if i!=j);H[j,j]=-sum(H[i,j] for i in range(n) if i!=j)
def predicates(M):
 out=[]
 for i,j,l in itertools.combinations(range(n),3):
  A,B,C=M[i,l],M[j,l],M[i,j]
  out.extend([A*(1+q[l]/(2*q[i]))-B*(1+q[l]/(2*q[j]))-C*(1-q[j]/q[i]),A*q[j]/q[i]-B])
 return out
orig=predicates(M);hs=predicates(H);repair=max([s.Rational(0)]+[a/(-b) for a,b in zip(orig,hs) if a>0])+s.Rational(1,10**30)
fixed=M+repair*H;assert all(z<0 for z in predicates(fixed))
conv=lambda z:mp.mpf(int(z.p))/int(z.q)
qm=list(map(conv,q))
def audit(M):
 S=mp.matrix([[conv(M[i,j])*mp.sqrt(qm[j]/qm[i]) for j in range(n)] for i in range(n)]);eig,U=mp.eigsy(S);w=[U[i,n-2]/mp.sqrt(qm[i]) for i in range(n)];Q=mp.mpf(0);coords=[]
 for j in range(n-1):
  Q+=qm[j];norm=mp.sqrt(Q*(Q+qm[j+1])/qm[j+1]);coords.append((sum(qm[i]*w[i] for i in range(j+1))-Q*w[j+1])/norm)
 if sum(coords)<0:w=[-z for z in w];coords=[-z for z in coords]
 return dict(arithmetic_margin=mp.nstr(w[0]+w[1]-2*w[2],80),gap=mp.nstr(eig[n-2]-eig[n-3],80),eigenvalues=[mp.nstr(z,50) for z in eig],perron=[mp.nstr(z,60) for z in coords],w=[mp.nstr(z,60) for z in w])
out=dict(original_search=p,original_exact_max_FG=str(max(orig)),repair=str(repair),unrepaired=audit(M),repaired=audit(fixed))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
