"""Independent direct-exponential and physical-parameter sensitivity checks."""
import itertools
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from activated_intermediate_model import ActivatedEvolution
from check_physical_rate_model import perturbed

ROOT=Path(__file__).resolve().parent


def independent(parameters,eta):
    mp.mp.dps=70
    x=[mp.mpf(str(z)) for z in parameters]
    E=[mp.mpf(0),x[0],x[0]+x[1],sum(x[:3])]
    pairs=list(itertools.combinations(range(4),2))
    B=[]
    for i,j in pairs:
        h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        B.append(E[j]+h)
    energies=E+B
    eta=mp.mpf(str(eta));weights=[mp.mpf(1)]*4+[eta]*6
    M=mp.matrix(10)
    for s,(i,j) in enumerate(pairs,4):
        for k in [i,j]:
            M[s,k]=2*mp.exp(-(energies[s]-energies[k]))
            M[k,s]=2/eta
    for j in range(10):M[j,j]=-sum(M[i,j] for i in range(10) if i!=j)
    def thermal(T):
        q=mp.matrix([g*mp.exp(-e/T) for g,e in zip(weights,energies)])
        return q/sum(q)
    pi=thermal(1);initial=[thermal(T) for T in x[-2:]]
    def distance(t,j):
        p=mp.expm(M*t)*initial[j]
        return sum(p[i]*mp.log(p[i]/pi[i]) for i in range(10))
    fast=ActivatedEvolution(parameters,float(eta))
    comparisons=[]
    for t in [0,30,60,120,220]:
        a=fast.distances(t)[:,0]
        b=[distance(t,j) for j in range(2)]
        err=max(abs(b[j]-mp.mpf(float(a[j]))) for j in range(2))
        comparisons.append(dict(time=t,KL=list(map(str,b)),double_error=str(err)))
    hits=[]
    for d in [mp.mpf('.001'),mp.mpf('.01')]:
        times=[]
        for j in range(2):
            a=mp.mpf(str(fast.hit(float(d),j)))
            # Independent secant iteration from the numerical estimate.
            r=mp.findroot(lambda t:distance(t,j)-d,(a-mp.mpf('.01'),a+mp.mpf('.01')),
                          tol=mp.mpf('1e-55'))
            times.append(r)
        hits.append(dict(threshold=str(d),times=list(map(str,times)),saving=str(1-times[1]/times[0])))
    # Exact block algebra separately proves the trace identity; numerical residual is a check.
    residual=max(abs(M[i,j]*pi[j]-M[j,i]*pi[i]) for i in range(10) for j in range(10))
    return dict(precision_digits=70,method="direct full ten-state matrix exponential",
                detailed_balance_residual=str(residual),time_checks=comparisons,hitting_checks=hits,
                interval_certificate=False)


if __name__=='__main__':
    x=json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']
    result=dict(nominal_eta=.05,independent=independent(x,.05))
    corners=list(itertools.product([-1,1],repeat=8))
    rng=np.random.default_rng(12102041)
    sample_records=[]
    for eta in [.025,.05,.1]:
        worst=1.;worst_x=None;fail=0;minimum_window=5.
        for signs in corners+list(rng.uniform(-1,1,(500,8))):
            y=perturbed(x,signs)
            m=ActivatedEvolution(y,eta).metrics()
            gain=m['minimum_time_saving_sampled']
            if gain<worst:worst=gain;worst_x=y.tolist()
            window=m['distance_window_tau_sampled']
            width=0 if window is None else window[1]-window[0]
            minimum_window=min(minimum_window,width)
            fail+=int(gain<.05 or width<.25)
        sample_records.append(dict(eta=eta,corners=256,random_points=500,
                                  minimum_time_saving_sampled=worst,
                                  minimum_distance_window_tau_sampled=minimum_window,
                                  failing_cases=fail,worst_time_parameters=worst_x))
        print(sample_records[-1],flush=True)
    result.update(robustness=sample_records,seed=12102041,
                  samples_are_not_a_certificate=True,
                  physical_scope="activated-intermediate mesoscopic extension, not continuum validation")
    (ROOT/'activated_intermediate_audit.json').write_text(json.dumps(result,indent=2)+'\n')
