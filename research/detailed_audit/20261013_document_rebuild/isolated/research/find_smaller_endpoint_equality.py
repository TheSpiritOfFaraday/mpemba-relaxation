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
 rad=min([p@p for p in cand if np.max(c+v@p)<1e-9]+[np.inf]);return rad

from scipy.optimize import minimize_scalar
cases=[]
for den in [2,3,4,5,6,8,10,20]:
 for m in range(2,12):
  for n in range(m+1,41):cases.append((n*np.log10(den),den,m,n))
for score,den,m,n in sorted(cases):
 if score>40*np.log10(4):break
 base=s.Rational(den-1,den);q0=np.array([1.,float(base),float(base)**m,float(base)**n]);vv=np.array([0,1,m,n],float);vv-=vv.mean();mat=np.array([q0[2:],vv[2:]])
 def objective(w2):
  ww=np.r_[1,w2,np.linalg.solve(mat,[-q0[0]-q0[1]*w2,-vv[0]-vv[1]*w2])]
  if ww[0]+ww[1]<2*ww[2] or ww[3]>0:return 1e4
  return min(1e4,feasible(q0,ww))
 opt=minimize_scalar(objective,bounds=(-1.1,0),method='bounded',options={'xatol':1e-6})
 if opt.fun<1-1e-6:
  w2=s.Rational(float(opt.x)).limit_denominator(10000);q=s.Matrix([1,base,base**m,base**n]);v=s.Matrix([0,1,m,n])-s.ones(4,1)*s.Rational(1+m+n,4)
  ww=s.Matrix([[q[2],q[3]],[v[2],v[3]]]).inv()*s.Matrix([-q[0]-q[1]*w2,-v[0]-v[1]*w2]);w=s.Matrix([1,w2,*ww]);out=exact_realize(q,w)
  print('FOUND',den,m,n,str(w2),opt.fun,out['feasible'],flush=True)
  if out['feasible']:
   final={'bath_weights':list(map(str,q)),'downward_rates':out['downward_rates'],'mode':list(map(str,w)),'energy_exponents':[0,1,m,n],'base':str(base),'exact_endpoint_covariance':'0','realization':out};path=Path(__file__).with_name('smaller_endpoint_equality.json');path.write_text(json.dumps(final,indent=2)+'\n');print('SAVED',len(path.read_text()),flush=True);break
