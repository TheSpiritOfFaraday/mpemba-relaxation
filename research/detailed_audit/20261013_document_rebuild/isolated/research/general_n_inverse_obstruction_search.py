"""Numerical PSD-dual/primal adversarial search for prescribed violating modes."""
import itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog,minimize
from scipy.linalg import null_space
rng=np.random.default_rng(123397);start=time.time();counts={'modes':0,'trace_dual':0,'dual_failure':0,'primal_attempts':0};best=-1
for trial in range(5000):
 n=5+trial%2;E=np.r_[0,np.cumsum(10**rng.uniform(-2,1,n-1))];q=np.exp(-E);q/=sum(q);Q2=sum(q[:2]);Q3=sum(q[:3]);threshold=(q[0]-q[1])*q[2]/(2*q[1]*Q3)
 coeff=np.r_[1,threshold*rng.uniform(.001,.999),10**rng.uniform(-5,5,n-3)];H=np.zeros((n,n-1))
 for j in range(n-1):H[:j+1,j]=1;H[j+1,j]=-sum(q[:j+1])/q[j+1]
 w=H@coeff;z=np.sqrt(q)*w;z/=np.linalg.norm(z);U=null_space(np.array([np.sqrt(q),z]));d=n-2;bpairs=list(itertools.combinations_with_replacement(range(d),2));Bs=[]
 for i,j in bpairs:
  B=np.zeros((d,d));B[i,j]=B[j,i]=1;Bs.append(B)
 Ms=[-np.sqrt(q[:,None]/q[None,:])*(U@B@U.T) for B in Bs];rows=[]
 for i,j,k in itertools.combinations(range(n),3):
  for typ in [0,1]:
   row=np.array([M[i,k]*(1+q[k]/(2*q[i]))-M[j,k]*(1+q[k]/(2*q[j]))-M[i,j]*(1-q[j]/q[i]) if typ==0 else M[i,k]*q[j]/q[i]-M[j,k] for M in Ms]);norm=max(abs(row));rows.append(row/norm)
 rows=np.array(rows);target=np.array([int(i==j) for i,j in bpairs]);sol=linprog(np.ones(len(rows)),A_eq=rows.T,b_eq=target,bounds=[(0,None)]*len(rows),method='highs');counts['modes']+=1
 if sol.success and np.max(abs(rows.T@sol.x-target))<1e-6:counts['trace_dual']+=1;continue
 counts['dual_failure']+=1
 if counts['primal_attempts']>=100:continue
 counts['primal_attempts']+=1
 def mat(x):return sum(B*x[j] for j,B in enumerate(Bs))
 x0=np.r_[target/d,-.1]
 constraints=[{'type':'eq','fun':lambda x:target@x[:-1]-1},{'type':'ineq','fun':lambda x:-rows@x[:-1]},{'type':'ineq','fun':lambda x:np.linalg.eigvalsh(mat(x[:-1]))-x[-1]}]
 op=minimize(lambda x:-x[-1],x0,method='SLSQP',constraints=constraints,bounds=[(-2,2)]*len(Bs)+[(-2,1)],options={'maxiter':200,'ftol':1e-11})
 mineig=min(np.linalg.eigvalsh(mat(op.x[:-1])));maxF=max(rows@op.x[:-1]);best=max(best,mineig)
 if mineig>1e-7 and maxF<1e-8:
  out=dict(candidate=True,counts=counts,q=q.tolist(),w=w.tolist(),coeff=coeff.tolist(),B=mat(op.x[:-1]).tolist(),U=U.tolist(),min_eig=float(mineig),max_constraint=float(maxF),elapsed=time.time()-start)
  Path('research/general_n_inverse_candidate.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));raise SystemExit
out=dict(candidate=False,counts=counts,best_primal_min_eigenvalue=float(best),elapsed=time.time()-start,seed=123397)
Path('research/general_n_inverse_obstruction_search_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
