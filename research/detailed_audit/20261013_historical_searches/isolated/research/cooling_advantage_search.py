"""Optimize cooling advantage using convex inverse-mode feasibility.
Discovery only. The full triplet predicate and eigenvalue assignment must be
checked independently before any proposed witness is called certified.
"""
import argparse,json,time,os
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution,brentq
from mode_realizability import realize,contrast_basis
from four_state_classifier import classify_mode,coordinates
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=3600);args=parser.parse_args()
P=Path(__file__).resolve().parent;start=time.time();calls=0;best=0.;bestdata=None;last=start
baseq=np.array([1,199/200,2/5,1/400]);baseE=-np.log(baseq);basew=np.array([1,-.75893729814493,-.61211930644096,-.00386630776426]);p=baseq/sum(baseq)
cc=contrast_basis(baseq).T@(p*basew);x0=np.r_[np.log(np.diff(baseE)),np.log(cc[:2]/cc[2])]
def evaluate(x,save=False):
 global calls,best,bestdata,last
 calls+=1;E=np.r_[0,np.cumsum(np.exp(x[:3]))]
 if E[-1]>20:return 10+E[-1]-20
 q=np.exp(-E);w=contrast_basis(q)@np.r_[np.exp(x[3:]),1.];w/=max(abs(w))
 r=realize(q,w,slack=1e-7)
 if not r['feasible']:
    rad=r['min_radius_squared'];return 1+min(100,max(0,rad-1))
 cl=classify_mode(E,w)
 if cl['cooling'] is not True:return max(0,min(cl['cooling_values']))+1e-5
 A,B,C=w[:3]-w[3]
 def phi(beta):
    u,v=coordinates(E,beta);return A+B*u+C*v
 if phi(0)<0:
    hot=0.;cold=brentq(phi,0.,1.,xtol=1e-13)
 else:
    mid=cl['vertex_beta']
    if mid is None or not 0<mid<1:return 1e-5
    hot=brentq(phi,0.,mid,xtol=1e-13);cold=brentq(phi,mid,1.,xtol=1e-13)
 def b(beta):
    zz=np.exp(-beta*E);return w@zz/sum(zz)
 bh,bc=b(hot),b(cold)
 if not (bc<0 and bh<0):return 1e-5
 gain=1-abs(bh)/abs(bc)
 if save or gain>best+1e-6:
    best=max(best,gain);bestdata=dict(E=E.tolist(),q=q.tolist(),w=w.tolist(),rates=r['rates'],
        hot_beta_limit=hot,cold_beta=cold,hot_overlap=bh,cold_overlap=bc,relative_advantage=gain,
        feasibility_radius_squared=r['min_radius_squared'],coordinates=x.tolist(),evaluations=calls,
        elapsed_seconds=time.time()-start,label='NUMERICAL DISCOVERY; beta=0 endpoint may require finite replacement')
    (P/'cooling_advantage_candidate.json').write_text(json.dumps(bestdata,indent=2)+'\n')
    print('BEST',gain,'beta',hot,cold,'E',E.tolist(),flush=True)
 return -gain

for run in range(10000):
 if time.time()-start>args.seconds:break
 # Alternate local energy neighborhoods and broad searches.
 if run%3==0:
    fixed=bestdata['coordinates'][:3] if bestdata else x0[:3]
    def fun(z):return evaluate(np.r_[fixed,z])
    result=differential_evolution(fun,[(-14,14)]*2,seed=111034+run,popsize=15,maxiter=500,tol=1e-8,polish=False)
 else:
    initial=np.array(bestdata['coordinates']) if bestdata else x0
    result=differential_evolution(evaluate,[(-8,2.5)]*3+[(-14,14)]*2,seed=111034+run,popsize=12,maxiter=800,tol=1e-7,polish=False,x0=initial)
 status=dict(pid=os.getpid(),status='running',run=run,elapsed_seconds=time.time()-start,evaluations=calls,best_advantage=best)
 (P/'cooling_advantage_status.json').write_text(json.dumps(status,indent=2)+'\n');print('RUN',run,json.dumps(status),flush=True)
status=dict(pid=os.getpid(),status='finished',elapsed_seconds=time.time()-start,evaluations=calls,best_advantage=best)
(P/'cooling_advantage_status.json').write_text(json.dumps(status,indent=2)+'\n');print('FINISHED',json.dumps(status),flush=True)
