"""Fixed-seed numerical checks, not a proof, of two finite-time entropy results."""
import json,itertools,time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm
from general_n_arithmetic_search import setup

SEED=202610112253
rng=np.random.default_rng(SEED)
results={'seed':SEED,'general_G':{},'tridiagonal':{}}

def record(stats,key,val,scale,context):
 s=stats.setdefault(key,dict(evaluations=0,accepted=0,skipped_tiny=0,negative=0,min_relative=1.0))
 s['evaluations']+=1
 if scale<1e-24:
  s['skipped_tiny']+=1;return
 s['accepted']+=1
 rel=float(val/scale)
 if rel<s['min_relative']:s['min_relative']=rel;s['worst']=context
 if val < -1e-10*scale:s['negative']+=1

def thermal(q,E,beta,side):
 p=np.exp(-beta*E-np.max(-beta*E));p/=sum(p)
 r=p/q;dr=side*r*(p@E-E)
 return r,dr

def run(kind,N):
 stats=results[kind];stats['generators']=N;stats['dimensions']={};stats['ratio_order_failure_count']=0
 stats['minimum_generator_offdiag']=float('inf');stats['largest_column_sum_residual']=0
 stats['expm_spotchecks']=[];stats['generators_also_satisfying_G']=0
 for case in range(N):
  n=3+case%6;stats['dimensions'][n]=stats['dimensions'].get(n,0)+1
  E=np.r_[0,np.cumsum(10**rng.uniform(-2,1,n-1))];E*=rng.uniform(1,15)/E[-1]
  q=np.exp(-E);q/=sum(q);sqrtq=np.sqrt(q)
  pairs,rows,Z,edges,d=setup(q)
  if kind=='general_G':
   K=np.zeros((n,n))
   for j in range(1,n):
    K[:j,j]=np.sort(10**rng.uniform(-3,3,j));K[j,:j]=K[:j,j]
   M=q[:,None]*K;np.fill_diagonal(M,-np.sum(M,axis=0))
   A=-np.tensordot(np.array([K[i,j] for i,j in pairs]),edges,axes=(0,0))
   assert np.min(A-np.diag(np.diag(A))) > -1e-9
   assert max(K[i,k]-K[j,k] for i,j,k in itertools.combinations(range(n),3))<=0
   zs=[1e-4,.05,.25,.5,.8,2/np.sqrt(5)]
   sides=[-1,1]
  else:
   off=10**rng.uniform(-2,2,n-2)
   A=np.diag(-10**rng.uniform(-2,2,n-1))+np.diag(off,1)+np.diag(off,-1)
   S=Z@A@Z.T;M=sqrtq[:,None]*S/sqrtq[None,:]
   reset=max(0,max(-M[i,j]/q[i] for i,j in pairs),max(-M[j,i]/q[j] for i,j in pairs))+1
   M+=reset*(q[:,None]*np.ones((1,n))-np.eye(n));A-=reset*np.eye(n-1)
   zs=[.001,.05,.3,1,3,8]
   sides=[1]
  stats['minimum_generator_offdiag']=min(stats['minimum_generator_offdiag'],float(min(M[i,j] for i in range(n) for j in range(n) if i!=j)))
  stats['largest_column_sum_residual']=max(stats['largest_column_sum_residual'],float(np.max(abs(M.sum(axis=0)))))
  assert np.max(abs(M*q[None,:]-(M*q[None,:]).T))<1e-9
  eigen,U=np.linalg.eigh(A);assert eigen[-1]<0
  rawgap=-eigen[-1]
  if max(M[i,k]*q[j]/q[i]-M[j,k] for i,j,k in itertools.combinations(range(n),3)) <=1e-10:stats['generators_also_satisfying_G']+=1
  A/=rawgap;eigen/=rawgap
  for z in zs:
   for side in sides:
    beta=1+side*z/E[-1]
    assert beta>0
    r0,dr0=thermal(q,E,beta,side)
    x0=Z.T@(sqrtq*(r0-1));y0=Z.T@(sqrtq*dr0)
    for tau in [0,.001,.01,.1,1,3,10]:
     exp=np.exp(eigen*tau);x=U@(exp*(U.T@x0));y=U@(exp*(U.T@y0))
     dev=(Z@x)/sqrtq;r=1+dev;dp=sqrtq*(Z@y)
     if np.min(r)<=0:raise AssertionError(('nonpositive evolved density',kind,case,z,tau,float(min(r))))
     ctx=dict(case=case,n=n,z=z,side=side,tau=tau,beta=beta,physical_time=tau/rawgap)
     derivs={'KL':np.log1p(dev)}
     if case%20==0 and z==zs[-1] and tau==.1 and side==1:
      P=expm(M*(tau/rawgap));pp=P@(q*r0);dd=P@(q*dr0)
      direct=float(dd@np.log(pp/q));spectral=float(dp@np.log1p(dev))
      h=1e-5/max(1,E[-1]);ds=[]
      for bb in [beta-h,beta+h]:
       rr,_=thermal(q,E,bb,1);pt=P@(q*rr);ds.append(float(pt@np.log(pt/q)))
      finite_difference=(ds[1]-ds[0])/(2*h)
      err=abs(direct-spectral)/max(abs(direct),1e-20)
      fderr=abs(finite_difference-direct)/max(abs(direct),1e-20)
      assert err<1e-7 and fderr<1e-5,(err,fderr)
      stats['expm_spotchecks'].append(dict(case=case,relative_derivative_error=err,finite_difference_relative_error=fderr,max_probability_error=float(max(abs(pp-q*r)))))
     if kind=='tridiagonal':
      derivs.update(chi2=2*dev,reverse_KL=dev/r,Hellinger=dev/(np.sqrt(r)*(np.sqrt(r)+1)),quartic=4*dev**3,total_variation=np.sign(dev))
     for name,grad in derivs.items():
      terms=dp*grad;record(stats,name,float(sum(terms)),float(sum(abs(terms))),ctx)
     if side==1 and min(x)>1e-12 and min(y)>1e-12:
      ratios=y/x
      if np.max(np.diff(ratios))>1e-8*max(abs(ratios)):
       stats['ratio_order_failure_count']+=1
       if 'first_ratio_order_failure' not in stats:
        stats['first_ratio_order_failure']=dict(**ctx,q=q.tolist(),E=E.tolist(),M=M.tolist(),A=A.tolist(),ratios=ratios.tolist())
  if (case+1)%50==0:print(kind,case+1,flush=True)

start=time.time();run('general_G',240);run('tridiagonal',240)
r=sp.Matrix([sp.Rational(3,5),sp.Rational(3,10),sp.Rational(1,10)]);z=sp.Matrix([1,2,-12]);D=sp.diag(*r)
J=r*sp.ones(1,3)-sp.eye(3)-D*z*z.T/(z.T*D*z)[0]
K21=sp.factor(J[0,1]-J[0,2]);assert K21==-sp.Rational(14,27)
assert J* r==sp.zeros(3,1);assert sp.ones(1,3)*J==sp.zeros(1,3)
assert J.eigenvals()=={sp.S(0):1,sp.S(-1):1,sp.S(-2):1}
A,B,C=J[0,2],J[1,2],J[0,1];t=r[1]/r[0];u=r[2]/r[1]
F=sp.factor(A*(1+t*u/2)-B*(1+u/2)-C*(1-t));G=sp.factor(A*t-B)
assert F==0 and G==-sp.Rational(2,9)
results['exact_J3']={'matrix':[[str(v) for v in J.row(i)] for i in range(3)],'K21':str(K21),'F123':str(F),'G123':str(G),'spectrum':['-2','-1','0']}
results['elapsed_seconds']=time.time()-start
results['tolerances']={'tiny_absolute_sum_cutoff':1e-24,'relative_negative_cutoff':1e-10,'ratio_order_cutoff':1e-8,'ratio_min_coordinates':1e-12}
results['ranges']={'n':[3,8],'energy_span':[1,15],'raw_general_reduced_rates':[.001,1000],'tridiagonal_raw_scales':[.01,100],'time_in_inverse_slowest_rate':[0,.001,.01,.1,1,3,10],'general_z':[.0001,.05,.25,.5,.8,2/np.sqrt(5)],'tridiagonal_heating_z':[.001,.05,.3,1,3,8]}
# A completed sampling run is not a passing check if it found a contradiction.
failures=[]
for kind in ('general_G','tridiagonal'):
 for metric,stats in results[kind].items():
  if isinstance(stats,dict) and stats.get('negative',0):
   failures.append(f'{kind}/{metric}: {stats["negative"]} negative derivatives')
if results['tridiagonal']['ratio_order_failure_count']:
 failures.append('Tridiagonal evolution violated the predicted contrast-ratio order')
results['passed']=not failures
results['failures']=failures
results['scope']='Fixed-seed floating-point sampling with declared tolerances and tiny-value exclusions; analytic proofs are separate.'
Path('research/finite_time_entropy_independent_audit.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
if failures: raise AssertionError('; '.join(failures))
