"""Standard-library exact certificate for a five-state triplet counterexample.
Adds E5=9000 log(10), beta_b=1, all downward rates from 5 equal to 1000.
Uses rational LDL inertia and an old-coordinate residual bound, avoiding a
spurious 1/sqrt(q5) loss. Thermal weights of state5 are exact rational powers.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json,time
from standalone_primary_check import Q,RATES,R,W_PREFIX,row_generator,thermal_root,require
start=time.time();n=5;A=F(1000);q=list(Q)+[F(1,10**9000)];old=row_generator(Q,RATES)
L=[[F(0) for _ in range(n)] for _ in range(n)]
for i in range(4):
 for j in range(4):L[i][j]=old[i][j]
 L[i][4]=A*q[4]/q[i];L[i][i]-=L[i][4];L[4][i]=A
L[4][4]=-4*A
require(all(sum(row)==0 for row in L),'row conservation')
require(all(L[i][j]>0 for i in range(n) for j in range(n) if i!=j),'strict positive rates')
require(all(q[i]*L[i][j]==q[j]*L[j][i] for i in range(n) for j in range(n)),'detailed balance')
for i,j,k in combinations(range(n),3):
 aa,bb,cc=L[k][i],L[k][j],L[j][i];t=q[j]/q[i];u=q[k]/q[j]
 FF=aa*(1+t*u/2)-bb*(1+u/2)-cc*(1-t);GG=aa*t-bb;kap=(aa+bb)/(cc*(1+t))
 require(FF<0 and GG<0,'all ten triplets strict')
 require(not ((kap>=1 and FF>0) or (kap<=1 and GG>0)),'original Boolean predicate')
def inertia(shift):
 mat=[[q[i]*(L[i][j]-(shift if i==j else 0)) for j in range(n)] for i in range(n)];positive=0
 for k in range(n):
  d=mat[k][k];require(d!=0,'nonzero exact LDL pivot');positive+=d>0
  for i in range(k+1,n):
   for j in range(k+1,n):mat[i][j]-=mat[i][k]*mat[k][j]/d
 return positive
counts=[]
for shift,count in [(-4001,5),(-3999,4),(-337,4),(-336,3),(-60,3),(-59,2),(-20,2),(-19,1)]:
 got=inertia(F(shift));require(got==count,'five-state exact spectral brackets');counts.append([shift,got])
w=list(W_PREFIX)+[-sum(q[i]*W_PREFIX[i] for i in range(3))/q[3]]
w.append(A*sum(w)/(R+4*A));mean=sum(qi*wi for qi,wi in zip(q,w))/sum(q);w=[wi-mean for wi in w]
require(sum(qi*wi for qi,wi in zip(q,w))==0,'five-state zero bath mean')
res=[sum(L[i][j]*w[j] for j in range(n))-R*w[i] for i in range(n)]
eta=F(1,10**7);gap=F(39)
require(sum(qi*ri*ri for qi,ri in zip(q,res))<eta**2,'weighted residual')
require(sum(qi*wi*wi for qi,wi in zip(q,w))>1,'trial weighted norm')
require(abs(res[4])<eta and abs(w[4])<1,'new-coordinate residual bound')
# The spectral theorem gives |lambda-R|<eta and old-coordinate error below
# 20 eta/gap. The new row equation supplies its own coordinate estimate.
olderr=20*eta/gap;newerr=(4*A*olderr+2*eta)/(4*A-20)
err=F(1,10**6);require(max(olderr,newerr)<err,'uniform coordinate error')
def overlap(root_n,new_root):
 roots=[thermal_root(qi,root_n) for qi in Q]+[(new_root,new_root)]
 lo=sum(wi*(a if wi>=0 else b) for wi,(a,b) in zip(w,roots));hi=sum(wi*(b if wi>=0 else a) for wi,(a,b) in zip(w,roots))
 zl=sum(a for a,b in roots);zh=sum(b for a,b in roots);corners=[a/b for a in [lo,hi] for b in [zl,zh]]
 return min(corners)-err,max(corners)+err
hot=overlap(1000,F(1,10**9));cold=overlap(100,F(1,10**90))
require(hot[1]<0 and cold[1]<0,'negative true overlaps')
margin=hot[0]-cold[1];require(margin>F(1,10000),'strict finite cooling advantage')
# Decimal output is explanatory only; all acceptance comparisons above exact.
out=dict(verified=True,new_energy='9000 log(10)',new_downward_rates='1000 to each original state',
 beta_hot='1/1000',beta_cold='1/100',beta_bath='1',triplets_checked=10,all_F_and_G_strict=True,
 spectral_counts=counts,uniform_overlap_error_bound=str(err),
 displayed_hot_interval=list(map(float,hot)),displayed_cold_interval=list(map(float,cold)),
 margin_proved_greater_than='1/10000',elapsed_seconds=time.time()-start,
 scope='Exact five-state extension of the original four-state witness; not a practical energy scale')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
