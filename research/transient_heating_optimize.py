"""Targeted differential-evolution search for finite-time KL heating.

The all-G conductance ordering is enforced by sorted logarithmic rates.
Scores are numerical discovery only; any negative score requires an independent
high-precision check.  No theorem is inferred from unsuccessful optimization.
"""
import argparse,json,time,os
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=3600)
args=parser.parse_args();P=Path(__file__).resolve().parent
start=time.time();calls=0;best=1.;bestdata=None;hits=0
def evaluate(x):
 global calls,best,bestdata,hits
 calls+=1
 E=np.r_[0,np.cumsum(np.exp(x[:3]))]
 if E[-1]>25:return 2+E[-1]-25
 q=np.exp(-E);q/=sum(q);K=np.zeros((4,4));offset=3
 for j in range(1,4):
  K[:j,j]=np.sort(np.exp(x[offset:offset+j]));K[j,:j]=K[:j,j];offset+=j
 L=K*q[None,:];np.fill_diagonal(L,-L.sum(axis=1))
 S=L*np.sqrt(q[:,None]/q[None,:]);val,V=np.linalg.eigh(S)
 if val[-2]>=-1e-8*max(1.,abs(val[0])):return 2
 beta=1+np.exp(x[9]);tt=np.exp(x[10])/(-val[-2])
 pp=np.exp(-beta*E);pp/=sum(pp);dpp=-pp*(E-pp@E)
 co=V[:,:3].T@((pp-q)/np.sqrt(q));dc=V[:,:3].T@(dpp/np.sqrt(q))
 decay=np.exp(tt*val[:3]);pt=q+np.sqrt(q)*(V[:,:3]@(co*decay))
 dt=np.sqrt(q)*(V[:,:3]@(dc*decay))
 if min(pt)<1e-12:return 2
 terms=dt*np.log(pt/q);scale=sum(abs(terms));deriv=sum(terms)
 if scale<1e-12:return 2
 score=deriv/scale
 if score<best-1e-7:
  best=score;bestdata=dict(E=E.tolist(),q=q.tolist(),K=K.tolist(),beta=float(beta),
   time=float(tt),scaled_time=float(np.exp(x[10])),score=float(score),derivative=float(deriv),
   p=pt.tolist(),eigenvalues=val.tolist(),coordinates=x.tolist(),evaluations=calls,
   elapsed_seconds=time.time()-start,label='numerical discovery; not a certificate')
  (P/'transient_heating_optimized_candidate.json').write_text(json.dumps(bestdata,indent=2)+'\n')
  print('BEST',score,'derivative',deriv,flush=True)
 if score<-1e-6:hits+=1
 return score
initial=None
old=P/'transient_heating_status.json'
if old.exists():
 b=json.loads(old.read_text())['best'];K=np.array(b['K'])
 initial=np.r_[np.log(np.diff(b['E'])),[np.log(K[i,j]) for j in range(1,4) for i in range(j)],
  np.log(b['beta']-1),np.log(b['scaled_time'])]
bounds=[(-7,2.5)]*3+[(-10,10)]*6+[(-9,5),(-10,1.5)]
for run in range(10000):
 if time.time()-start>args.seconds or hits:break
 if bestdata is not None:initial=bestdata['coordinates']
 result=differential_evolution(evaluate,bounds,seed=111100+run,popsize=12,maxiter=500,
  tol=1e-8,polish=False,x0=initial)
 status=dict(pid=os.getpid(),run=run,status='running',elapsed_seconds=time.time()-start,
  evaluations=calls,best_score=best,hits=hits)
 (P/'transient_heating_optimized_status.json').write_text(json.dumps(status,indent=2)+'\n')
 print('RUN',json.dumps(status),flush=True)
status.update(status='finished',elapsed_seconds=time.time()-start,evaluations=calls,best_score=best,hits=hits)
(P/'transient_heating_optimized_status.json').write_text(json.dumps(status,indent=2)+'\n')
print('FINISHED',json.dumps(status),flush=True)
