"""Outward interval certificate for a fixed decimal Arrhenius candidate.

Uniformization: exp(Mt)p = exp(-rho*t) sum_k (rho*t)^k/k! P^k p.
P is stochastic for the exact model. Every exact component of P^k p lies
in [0,1], so the remaining Poisson mass bounds each component's tail.
Monotonicity of relative entropy under the Markov semigroup then covers
whole time cells and whole threshold cells, not just their endpoints.
mpmath.iv supplies outward-rounded elementary interval operations.
This is an executable interval certificate, not a Lean formalization.
"""
import argparse
import itertools
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from physical_rate_model import Evolution

mp.iv.dps=60
I=mp.iv.mpf


def lower(a):return float(np.nextafter(float(a.a),-np.inf))
def upper(a):return float(np.nextafter(float(a.b),np.inf))
def strictly_positive(a):return bool(a.a>0)


def physical_box(parameters):
    """Certify the stated parameter limits for the entire +/-5% box."""
    x=[I(str(v))*I(['.95','1.05']) for v in parameters[:6]]
    temps=[1+(I(str(v))-1)*I(['.95','1.05']) for v in parameters[-2:]]
    assert all(bool(g>=I('.1')) for g in x[:3])
    assert bool(sum(x[:3])<=5)
    assert all(bool(h>=2) and bool(h<=6) for h in x[3:])
    assert bool(temps[0]>=I('1.25')) and bool(temps[1]<=10)
    assert bool(temps[1]>temps[0])
    rates=[]
    for i,j in itertools.combinations(range(4),2):
        h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        gap=sum(x[i:j])
        rates.extend([mp.iv.exp(-h),mp.iv.exp(-h-gap)])
    rate_min=min(a.a for a in rates);rate_max=max(a.b for a in rates)
    ratio=rate_max/rate_min
    assert bool(ratio<1000)
    return dict(scope="entire independent +/-5% parameter box",
                energy_span_upper=upper(sum(x[:3])),
                min_adjacent_gap_lower=min(lower(g) for g in x[:3]),
                barrier_excess_range=[min(lower(h) for h in x[3:]),max(upper(h) for h in x[3:])],
                rate_ratio_upper=upper(ratio),
                preparation_ranges=[[lower(t),upper(t)] for t in temps],
                performance_certified=False)


class Certificate:
    def __init__(self,parameters):
        x=[I(str(v)) for v in parameters]
        self.E=[I(0),x[0],x[0]+x[1],sum(x[:3])]
        self.M=[[I(0) for _ in range(4)] for _ in range(4)]
        for i,j in itertools.combinations(range(4),2):
            h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
            self.M[i][j]=mp.iv.exp(-h)
            self.M[j][i]=mp.iv.exp(-h-(self.E[j]-self.E[i]))
        for j in range(4):
            self.M[j][j]=-sum(self.M[i][j] for i in range(4) if i!=j)
        self.rho=I('.5')
        assert all(bool(-self.M[j][j]<self.rho) for j in range(4))
        self.P=[[I(int(i==j))+self.M[i][j]/self.rho for j in range(4)] for i in range(4)]
        assert all(bool(a>=0) for row in self.P for a in row)
        self.pi=self.thermal(I(1))
        self.initial=[self.thermal(x[-2]),self.thermal(x[-1])]
        self.max_terms=0
        self.max_tail=0.

    def thermal(self,T):
        q=[mp.iv.exp(-e/T) for e in self.E]
        z=sum(q)
        return [a/z for a in q]

    def inertia_at(self,shift):
        q=[mp.iv.exp(-e) for e in self.E]
        S=[[self.M[i][j]*mp.iv.sqrt(q[j]/q[i])-I(str(shift))*int(i==j)
            for j in range(4)] for i in range(4)]
        negative=0
        for k in range(4):
            pivot=S[k][k]
            assert bool(pivot<0) or bool(pivot>0),('uncertain LDL pivot',k)
            negative+=int(bool(pivot<0))
            for i in range(k+1,4):
                for j in range(k+1,4):
                    S[i][j]-=S[i][k]*S[k][j]/pivot
        return negative

    def distances(self,time):
        t=I(str(time));mu=self.rho*t
        v=[row[:] for row in self.initial]
        weight=mp.iv.exp(-mu)
        accum=[[weight*a for a in row] for row in v]
        k=0
        while True:
            next_weight=weight*mu/(k+1)
            if bool(mu<I(k+2)):
                tail=next_weight/(1-mu/(k+2))
                if bool(tail<I('1e-45')):break
            k+=1
            assert k<3000
            v=[[sum(self.P[i][j]*row[j] for j in range(4)) for i in range(4)] for row in v]
            weight=next_weight
            accum=[[a+weight*b for a,b in zip(row,new)] for row,new in zip(accum,v)]
        self.max_terms=max(self.max_terms,k+1)
        self.max_tail=max(self.max_tail,upper(tail))
        remainder=I([0,tail.b])
        vals=[]
        for row in accum:
            p=[a+remainder for a in row]
            assert all(bool(a>0) for a in p)
            vals.append(sum(a*mp.iv.log(a/b) for a,b in zip(p,self.pi)))
        return vals


def run(parameters):
    cert=Certificate(parameters)
    physical_box_bounds=physical_box(parameters)
    evo=Evolution(parameters)
    left=-1.1/evo.tau;right=-.9/evo.tau
    assert cert.inertia_at(left)==2
    assert cert.inertia_at(right)==3
    tau_lo=1/I(str(-left));tau_hi=1/I(str(-right))
    # The reported fixed interval is chosen inside the wide sampled reversal.
    # Time units are fixed by the nominal double tau; certification uses
    # the actual decimal times printed here, not an unverified eigenvalue.
    t1=.4*evo.tau;t2=.8*evo.tau
    assert strictly_positive(I(str(t2))-I(str(t1))-I('.25')*tau_hi)
    assert strictly_positive(5*tau_lo-I(str(t2)))
    times=np.linspace(t1,t2,81)
    D=[cert.distances(float(t)) for t in times]
    margins=[]
    for i in range(len(times)-1):
        margin=I('.9')*D[i+1][0]-D[i][1]
        assert strictly_positive(margin),('distance cell',i,str(margin))
        margins.append(lower(margin))
    assert strictly_positive(D[-1][0]-I('.001'))
    thresholds=np.geomspace(.001,.01,61)
    roots=[]
    for d in thresholds:
        brackets=[]
        for j in range(2):
            approx=evo.hitting(float(d),j)*evo.tau
            lo=float(approx-1e-7);hi=float(approx+1e-7)
            assert strictly_positive(cert.distances(lo)[j]-I(str(float(d))))
            assert strictly_positive(I(str(float(d)))-cert.distances(hi)[j])
            brackets.append([lo,hi])
        roots.append(brackets)
    saving_bounds=[]
    for i in range(len(thresholds)-1):
        # hot time is largest at the lower threshold, cold time is
        # smallest at the upper threshold.
        bound=1-I(str(roots[i][1][1]))/I(str(roots[i+1][0][0]))
        assert strictly_positive(bound-I('.05')),('threshold cell',i,str(bound))
        saving_bounds.append(lower(bound))
    initial=cert.distances(0)
    assert strictly_positive(initial[1]-initial[0])
    assert strictly_positive(initial[0]-I('.01'))
    return dict(parameters=parameters,arithmetic="mpmath.iv, 60 decimal digits",
                uniformization_rate='.5',max_terms=cert.max_terms,max_tail=cert.max_tail,
                time_interval_attempt_units=[t1,t2],
                time_interval_nominal_tau_units=[.4,.8],
                slow_eigenvalue_bracket=[left,right],
                tau_interval_attempt_units=[lower(tau_lo),upper(tau_hi)],
                distance_cell_count=80,minimum_KL_inequality_margin=min(margins),
                cold_distance_floor_lower=lower(D[-1][0]),
                threshold_band=[.001,.01],threshold_cells=60,
                certified_min_time_saving_lower=min(saving_bounds),
                threshold_root_brackets=roots,thresholds=thresholds.tolist(),
                parameter_neighborhood_certified=False,
                parameter_box_physical_constraints=physical_box_bounds,
                scope="fixed decimal model; whole time interval and whole threshold band")


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('receipt')
    args=parser.parse_args()
    source=Path(args.receipt)
    x=json.loads(source.read_text())['result']['parameters']
    result=run(x)
    source.with_name(source.stem+'_interval_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['threshold_root_brackets','thresholds']},indent=2))
