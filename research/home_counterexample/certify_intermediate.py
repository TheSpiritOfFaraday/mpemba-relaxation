"""Finite 8-percent strict counterexample, with exact arithmetic acceptance.
Reuses the archived symbolic family expressions, independently verifies their
specialization by exact matrix identities and inertia. Thermal intervals use
integer-power comparisons; perturbation bound follows the stated spectral lemma.
No claim of experimental robustness or useful finite-time signal.
"""
import sys,json
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from standalone_primary_check import approximate_overlap,positive_inertia,require
from check_intermediate import parse,source
import itertools
q=[F(1),F(2999,3000),F(3,5),F(1,10**6)];e=F(1,1000)
w=[F(1),F(-4,5),-F(1,3)-4*e/9+5*e**3/3,-e]
pi=[x/sum(q) for x in q]
M=[[parse(s,e) for s in row] for row in source['generator']]
require(sum(pi[i]*w[i] for i in range(4))==0,'bath orthogonality')
require(all(sum(w[i]*M[i][j] for i in range(4))==-3*w[j] for j in range(4)),'exact eigenmode -3')
require(all(sum(M[i][j] for i in range(4))==0 for j in range(4)),'conservation')
require(all(M[i][j]*pi[j]==M[j][i]*pi[i] for i in range(4) for j in range(4)),'reversibility')
require(all(M[i][j]>0 for i in range(4) for j in range(4) if i!=j),'positive rates')
L=[list(row) for row in zip(*M)]
count,pivots=positive_inertia(L,pi,F(-10,3))
require(count==2,'only zero and slow eigenvalue above -10/3; gap exceeds 1/3')
H=[[F(0) for _ in range(4)] for _ in range(4)]
for (i,j),k in zip(itertools.combinations(range(4),2),[100,10,1,11,2,3]):H[i][j]=pi[i]*k;H[j][i]=pi[j]*k
for j in range(4):H[j][j]=-sum(H[i][j] for i in range(4) if i!=j)
perturbed=[[M[i][j]+e**4*H[i][j] for j in range(4)] for i in range(4)]
res=[]
for i,j,k in itertools.combinations(range(4),3):
 A,B,C=perturbed[i][k],perturbed[j][k],perturbed[i][j];t=q[j]/q[i];u=q[k]/q[j]
 ff=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);gg=A*t-B;kappa=(A+B)/(C*(1+t))
 require(ff<0 and gg<0,'strict F/G')
 require(not ((kappa>=1 and ff>0) or (kappa<=1 and gg>0)),'original predicate fails')
 res.append(dict(triplet=[i+1,j+1,k+1],F=str(ff),G=str(gg)))
# H's reversible Dirichlet form has operator norm <= max reduced conductance=100.
alpha=100*e**4;gap=F(1,3)
require(gap-2*alpha>F(1,4),'perturbed slow remains isolated')
require(sum(pi[i]*w[i]**2 for i in range(4))<=1,'weighted norm <=1')
require(min(pi)>=e**2/3,'coordinate norm conversion')
# sin(angle)<=alpha/(gap-alpha); oriented unit-vector difference<=sqrt(2)*sin.
# sqrt(2)*sqrt(3)<3, hence sup error < 3*alpha/((gap-alpha)*e).
err=1000*e**3
require(3*alpha/((gap-alpha)*e)<err,'uniform overlap error bound')
hot=approximate_overlap(q,w,10000)
cold=approximate_overlap([x**9 for x in q],w,100)
require(hot[1]+err<0 and cold[1]+err<0,'negative perturbed overlaps')
lo=1-(-hot[0]+err)/(-cold[1]-err)
hi=1-(-hot[1]-err)/(-cold[0]+err)
require(lo>F(8,100),'strict amplitude advantage exceeds 8 percent')
out=dict(verified=True,epsilon=str(e),perturbation=str(e**4),q=list(map(str,q)),beta_hot='1/10000',beta_cold='9/100',beta_bath='1',generator=[[str(x) for x in r] for r in perturbed],triplet_residuals=res,unperturbed_eigenfunction=list(map(str,w)),unperturbed_LDL_pivots=list(map(str,pivots)),unperturbed_spectral_gap_lower='1/3',perturbed_spectral_gap_lower='1/4',uniform_overlap_error=str(err),amplitude_advantage_enclosure=[str(lo),str(hi)],display_percent=[100*float(lo),100*float(hi)],scope='Exact finite mathematical witness using spectral perturbation bound; no physical implementation or tolerance-box certification')
(Path(__file__).parent/'intermediate_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['verified','epsilon','perturbation','display_percent','uniform_overlap_error','scope']},indent=2))
