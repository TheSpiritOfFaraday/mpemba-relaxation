"""Exact bath-collision realization and independent event-clock simulation.

Each unordered pair has a Poisson contact clock of rate nu0=1. A fresh
bath energy Y~Exp(beta_b) is sampled. A system in endpoint j crosses to i
iff E_j+Y>=B_ij, after which Y'=Y+E_j-E_i. Rejected contacts do nothing.
Integrating Y gives exactly the report's Arrhenius rates.

Simulation thins the barrier-excess gate analytically, retaining its
waiting-time cost. Conditional on Y>=h, Y=h+Exp(beta_b), so only an
uphill residual-energy test remains. No fixed time step is used.
"""
import itertools
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm
from physical_rate_model import construct,thermal

ROOT=Path(__file__).resolve().parent

def algebra(parameters):
    beta,d,y,h,e=sp.symbols('beta d y h e',positive=True)
    acceptance=sp.integrate(beta*sp.exp(-beta*y),(y,h+d,sp.oo))
    assert sp.simplify(acceptance-sp.exp(-beta*(h+d)))==0
    assert sp.simplify(sp.integrate(beta*sp.exp(-beta*y),(y,h,sp.oo))-sp.exp(-beta*h))==0
    # Low state e, high state e+d; bath energy shifts by -d on the upward move.
    assert sp.simplify(sp.exp(-beta*(e+y))-sp.exp(-beta*(e+d+y-d)))==0
    x=[Fraction(str(v)) for v in parameters]
    E=[Fraction(0),x[0],x[0]+x[1],sum(x[:3])]
    # Binary labels: 00,10,11,01. Basins are fixed first coordinate.
    epsilon_x,epsilon_y,J=E[1],E[3],E[2]-E[1]-E[3]
    labels=[(0,0),(1,0),(1,1),(0,1)]
    assert [epsilon_x*a+epsilon_y*b+J*a*b for a,b in labels]==E
    checks=0
    for i,j in itertools.combinations(range(4),2):
        excess=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        barrier=E[j]+excess
        for state in [i,j]:
            threshold=barrier-E[state]
            for bath in [Fraction(0),threshold/2,threshold,threshold+Fraction(1,10),threshold*2]:
                def collision(s,Y):
                    if E[s]+Y<barrier:return s,Y
                    target=j if s==i else i
                    return target,Y+E[s]-E[target]
                target,outgoing=collision(state,bath)
                assert outgoing>=0
                assert E[state]+bath==E[target]+outgoing
                assert collision(target,outgoing)==(state,bath)
                checks+=1
    return dict(symbolic_acceptance=str(acceptance),exact_energy_and_involution_checks=checks,
                binary_Hamiltonian={'epsilon_x':str(epsilon_x),'epsilon_y':str(epsilon_y),'J':str(J)},
                theorem='For every allowed symmetric saddle and beta>0, a_ij=nu0 exp[-beta(B_ij-E_j)].',
                bath_assumptions=['independent fresh canonical energy on each contact','constant density of bath energy states','independent Poisson pair-contact clocks','instantaneous idealized energy-conserving contact map'])

def event_audit(parameters,N=500000,seed=12102043):
    E,B,M=construct(parameters);pairs=np.array(list(itertools.combinations(range(4),2)))
    excess=np.array([B[i,j]-E[j] for i,j in pairs])
    eligible_rates=np.exp(-excess);total=sum(eligible_rates)
    rng=np.random.default_rng(seed)
    pi=thermal(E,1);times=[0.,30.,60.,120.,210.]
    records=[];maxerror=0.;eligible_events=0;accepted_events=0
    for prep,T in enumerate(parameters[-2:]):
        initial=thermal(E,T)
        state=rng.choice(4,N,p=initial);last=0.
        for t in times:
            count=rng.poisson(total*(t-last),N)
            eligible_events+=int(count.sum())
            for step in range(int(count.max(initial=0))):
                ids=np.flatnonzero(count>step)
                pair=rng.choice(6,len(ids),p=eligible_rates/total)
                i,j=pairs[pair].T;s=state[ids]
                incident=(s==i)|(s==j)
                target=np.where(s==i,j,i)
                residual=rng.exponential(1.,len(ids))
                accept=incident&(residual>=E[target]-E[s])
                state[ids[accept]]=target[accept]
                accepted_events+=int(accept.sum())
            empirical=np.bincount(state,minlength=4)/N
            exact=expm(M*t)@initial
            error=float(np.max(np.abs(empirical-exact)));maxerror=max(maxerror,error)
            # Per-coordinate diagnostic envelope; no coverage theorem is inferred.
            bound=6*np.sqrt(exact*(1-exact)/N)+2/N
            assert np.all(np.abs(empirical-exact)<bound),(prep,t,empirical,exact)
            kl=float(sum(empirical*np.log(empirical/pi)))
            exactkl=float(sum(exact*np.log(exact/pi)))
            records.append(dict(preparation=['cold','hot'][prep],time=t,empirical_probability=empirical.tolist(),
                                direct_exponential_probability=exact.tolist(),maximum_probability_error=error,
                                empirical_KL=kl,exact_KL=exactkl,
                                empirical_mean_energy=float(empirical@E),exact_mean_energy=float(exact@E),
                                bath_mean_energy=float(pi@E)))
            last=t
    return dict(trajectories_per_preparation=N,seed=seed,clock='Poisson contacts including rejected-gate waiting time by exact thinning',
                fixed_time_step=False,eligible_contacts=eligible_events,accepted_transitions=accepted_events,
                maximum_probability_error=maxerror,all_diagnostic_checks_passed=True,records=records,
                simulation_is_not_a_robustness_certificate=True)

if __name__=='__main__':
    x=json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']
    result=dict(parameters=x,algebra=algebra(x),event_audit=event_audit(x),
                scope='Exact stochastic collision realization of the same four-state rates, not an autonomous Hamiltonian or identified material.')
    (ROOT/'thermal_collision_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='event_audit'},indent=2))
    print('Maximum probability error:',result['event_audit']['maximum_probability_error'])
