"""Independent 80-digit direct-exponential audit and sampled robustness.

The independent construction uses saddle energies directly, not the
symmetric eigendecomposition in physical_rate_model. No interval claims.
"""
import argparse
import itertools
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from physical_rate_model import Evolution

ROOT=Path(__file__).resolve().parent


def perturbed(x, signs, radius=.05):
    y=np.array(x)*(1+radius*np.array(signs))
    y[-2:]=1+(np.array(x[-2:])-1)*(1+radius*np.array(signs[-2:]))
    return y


def high_precision(x):
    mp.mp.dps=80
    x=list(map(lambda v:mp.mpf(str(v)),x))
    E=[mp.mpf(0),x[0],x[0]+x[1],sum(x[:3])]
    M=mp.matrix(4)
    for i,j in itertools.combinations(range(4),2):
        h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        saddle=max(E[i],E[j])+h
        M[i,j]=mp.exp(-(saddle-E[j]))
        M[j,i]=mp.exp(-(saddle-E[i]))
    for j in range(4): M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
    def thermal(T):
        q=mp.matrix([mp.exp(-e/T) for e in E])
        return q/sum(q)
    pi=thermal(1)
    initials=[thermal(T) for T in x[-2:]]
    def distance(t,j):
        p=mp.expm(M*t)*initials[j]
        return sum(p[i]*mp.log(p[i]/pi[i]) for i in range(4))
    tau=Evolution(list(map(float,x))).tau
    comparisons=[]
    maxerr=0
    for s in [0,.1,.25,.5,1,2]:
        t=mp.mpf(str(s))*mp.mpf(str(tau))
        vals=[distance(t,j) for j in range(2)]
        fast=Evolution(list(map(float,x))).distances(s)[:,0]
        maxerr=max(maxerr,max(abs(mp.mpf(float(fast[j]))-vals[j]) for j in range(2)))
        comparisons.append(dict(time_tau=s, KL_cold=str(vals[0]), KL_hot=str(vals[1])))
    thresholds=[]
    for d in [mp.mpf('.001'),mp.mpf('.0031622776601683793319988935444327185337'),mp.mpf('.01')]:
        times=[]
        for j in range(2):
            # Bisection relies on KL contractivity of this irreducible Markov semigroup.
            lo=mp.mpf(0);hi=mp.mpf(str(20*tau))
            for _ in range(160):
                mid=(lo+hi)/2
                if distance(mid,j)>d:lo=mid
                else:hi=mid
            times.append((lo+hi)/2)
        thresholds.append(dict(delta=str(d),times=list(map(str,times)),
                               saving=str(1-times[1]/times[0])))
    trips=[]
    q=[mp.exp(-e) for e in E]
    for i,j,k in itertools.combinations(range(4),3):
        A=M[i,k];B=M[j,k];C=M[i,j];t=q[j]/q[i];u=q[k]/q[j]
        F=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);G=A*t-B
        trips.append(dict(triplet=[i+1,j+1,k+1], F=str(F),G=str(G),
                          failure=bool(F<=0 and G<=0)))
    return dict(precision_digits=80,method="direct mpmath matrix exponential",
                detailed_balance_residual=str(max(abs(M[i,j]*pi[j]-M[j,i]*pi[i]) for i in range(4) for j in range(4))),
                max_double_vs_direct_KL_error=str(maxerr), comparisons=comparisons,
                hitting_checks=thresholds, triplets=trips,
                all_original_triplet_tests_fail=all(t['failure'] for t in trips))


def sample_box(x,random_count):
    rng=np.random.default_rng(12102027)
    directions=list(itertools.product([-1,1],repeat=8))
    directions.extend(rng.uniform(-1,1,(random_count,8)))
    worst_time=1.;max_rate=0.;worst_window=5.;worst_case=None
    envelope_low=np.ones(31)
    envelope_high=-np.ones(31)
    failures=0
    for direction in directions:
        y=perturbed(x,direction)
        m=Evolution(y).metrics()
        gain=m['min_time_saving_sampled']
        tt=np.array(m['hitting_times_tau'])
        savings=1-tt[1]/tt[0]
        envelope_low=np.minimum(envelope_low,savings)
        envelope_high=np.maximum(envelope_high,savings)
        window=m['distance_window_sampled']
        width=0 if window is None else window[1]-window[0]
        if gain<worst_time:worst_time=gain;worst_case=y.tolist()
        worst_window=min(worst_window,width)
        max_rate=max(max_rate,m['rate_ratio'])
        if gain<.05 or width<.25 or m['rate_ratio']>1000:failures+=1
    return dict(corners=256,random_interior=random_count,seed=12102027,
                min_time_saving_sampled=worst_time,min_window_width_tau_sampled=worst_window,
                max_rate_ratio_sampled=max_rate,failing_cases=failures,
                thresholds=m['thresholds'],time_saving_envelope_low=envelope_low.tolist(),
                time_saving_envelope_high=envelope_high.tolist(),
                worst_time_parameters=worst_case,certified_neighborhood=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('receipt')
    parser.add_argument('--samples',type=int,default=1000)
    args=parser.parse_args()
    source=Path(args.receipt)
    x=json.loads(source.read_text())['result']['parameters']
    result=dict(source=str(source),parameters=x,nominal=Evolution(x).metrics(),
                independent=high_precision(x),robustness=sample_box(x,args.samples))
    dest=source.with_name(source.stem+'_audit.json')
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['robustness'],indent=2))
    print('All triplet tests fail:',result['independent']['all_original_triplet_tests_fail'])
