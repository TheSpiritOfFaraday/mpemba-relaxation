import numpy as np,scipy.linalg as la,itertools,json,sys
import sympy as s
from exact_mode_realizability import exact_realize
from pathlib import Path
pairs=list(itertools.combinations(range(4),2));triples=list(itertools.combinations(range(4),3))
bs=[np.eye(2)/2,np.diag([1,-1])/2,np.array([[0,1],[1,0]])/2]
def feasible(q,w):
 q=q/q.sum();U=la.null_space(np.vstack([q,q*w]));down=np.array([[-q[i]*(U[i]@B@U[j]) for B in bs] for i,j in pairs]);rows=[]
 for i,j,k in triples:
  A,B,C=[down[pairs.index(p)] for p in [(i,k),(j,k),(i,j)]]
  rows.extend([A*(1+q[k]/(2*q[i]))-B*(1+q[k]/(2*q[j]))-C*(1-q[j]/q[i]),A*q[j]/q[i]-B])
 rows=np.array(rows);rows/=np.linalg.norm(rows,axis=1)[:,None];c=rows[:,0];v=rows[:,1:];cand=[np.zeros(2)]
 for cc,vv in zip(c,v):
  if vv@vv>1e-20:cand.append(-cc*vv/(vv@vv))
 for i,j in itertools.combinations(range(8),2):
  V=v[[i,j]]
  if abs(np.linalg.det(V))>1e-12:cand.append(np.linalg.solve(V,-c[[i,j]]))
 rad=min([p@p for p in cand if np.max(c+v@p)<1e-9]+[np.inf]);return rad<1-1e-7
candidates=[];tested=0
for den in [2,3,4,5,6,10,20]:
 base=s.Rational(den-1,den)
 for m in range(2,16):
  for n in range(m+1,61):
   q=np.array([1.,float(base),float(base)**m,float(base)**n]);v=np.array([0,1,m,n],float);v-=v.mean()
   for w2 in [s.Rational(-99,100),s.Rational(-19,20),s.Rational(-9,10),s.Rational(-4,5),s.Rational(-1,2),s.Rational(0)]:
    w=np.r_[1,float(w2),np.linalg.solve(np.array([q[2:],v[2:]]),[-q[0]-q[1]*float(w2),-v[0]-v[1]*float(w2)])];tested+=1
    if w[0]+w[1]<2*w[2]-1e-8 or w[3]>1e-8:continue
    if feasible(q,w):
     score=n*len(str(den));candidates.append((score,den,m,n,str(w2)))
print('Numeric candidates',len(candidates),'tested',tested,flush=True)
for score,den,m,n,w2 in sorted(candidates)[:30]:
 base=s.Rational(den-1,den);q=s.Matrix([1,base,base**m,base**n]);mm=s.Matrix([0,1,m,n]);v=mm-s.ones(4,1)*sum(mm)/4;w2=s.Rational(w2)
 ww=s.Matrix([[q[2],q[3]],[v[2],v[3]]]).inv()*s.Matrix([-q[0]-q[1]*w2,-v[0]-v[1]*w2]);w=s.Matrix([1,w2,*ww])
 result=exact_realize(q,w)
 print('Exact',den,m,n,w2,result['feasible'],flush=True)
 if result['feasible']:
  out={'bath_weights':list(map(str,q)),'downward_rates':result['downward_rates'],'mode':list(map(str,w)),'energy_exponents':[0,1,m,n],'base':str(base),'exact_endpoint_covariance':'0','realization':result}
  path=Path(__file__).with_name('small_endpoint_equality.json');path.write_text(json.dumps(out,indent=2)+'\n');print('SAVED',len(path.read_text()),'bytes',flush=True);break
