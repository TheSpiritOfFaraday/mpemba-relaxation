"""Adversarial numerical and exact-polynomial checks of the finite classifier."""
import json,time,itertools
from pathlib import Path
import numpy as np
import sympy as sy
import mpmath as mp
from four_state_classifier import classify_mode,coordinates,slope_x
from exact_polynomial_reference import sign_changes_open
ROOT=Path(__file__).resolve().parents[1];rng=np.random.default_rng(111028);start=time.time()
# Independent derivative using normalized probabilities and covariance, not u/v.
def raw_derivative(E,w,bs):
    # Sign-equivalent derivative numerator, factored by exp(-beta*E2).
    # Independent pairwise covariance formula prevents underflow at low T.
    return sum((w[i]-w[j])*(E[j]-E[i])*np.exp(-bs*(E[i]+E[j]-E[1]))
               for i in range(4) for j in range(i+1,4))
counts={'random':0,'strict_results':0,'boundary_or_uncertain':0,'grid_missed_classified_effect':0};examples={};errors=[]
for n in range(10000):
    E=np.r_[0,np.cumsum(10**rng.uniform(-3,1,3))];w=rng.normal(size=4)
    if n%2==0:
        q=np.exp(-E);w-=w@q/sum(q)
    out=classify_mode(E,w)
    for side,bs in [('cooling',np.linspace(0,1,801)),('heating',np.r_[np.linspace(1,5,401),np.geomspace(5, max(6,40/min(np.diff(E))),400)])]:
        pred=out[side];der=raw_derivative(E,w,bs)
        scale=max(abs(der));seen=(min(der)<-1e-9*scale and max(der)>1e-9*scale)
        counts['random']+=1
        if pred=='boundary_or_uncertain':counts['boundary_or_uncertain']+=1;continue
        counts['strict_results']+=1
        if seen and pred is False:errors.append({'E':E.tolist(),'w':w.tolist(),'side':side,'out':out,'issue':'classifier missed grid sign change'})
        if pred is True and not seen:
            counts['grid_missed_classified_effect']+=1
            # Evaluate independently at the theorem's unique interior extremum.
            if out['vertex_beta'] is not None and ((side=='cooling' and out['vertex_beta']<1) or (side=='heating' and out['vertex_beta']>1)):
                bs2=np.r_[bs,out['vertex_beta']];dd=raw_derivative(E,w,bs2)
                if not(min(dd)<0<max(dd)):errors.append({'E':E.tolist(),'w':w.tolist(),'side':side,'issue':'unconfirmed predicted sign change'})
            else:
                errors.append({'E':E.tolist(),'w':w.tolist(),'side':side,'issue':'unconfirmed endpoint prediction'})
    key=str((out['cooling'],out['heating']))
    if key not in examples:examples[key]={'E':E.tolist(),'w':w.tolist(),'classification':out}
    if n%1000==0:print('random progress',n,counts,'errors',len(errors),flush=True)

# Exact polynomial reference, equal energy gaps with q_i=(1/2)^i.
x=sy.Symbol('x');Z=1+x+x*x+x**3;polychecks=0;boundary=0
for n in range(300):
    ww=list(map(int,rng.integers(-20,21,4)))
    if len(set(ww))==1:continue
    N=sum(ww[i]*x**i for i in range(4));H=sy.Poly(sy.expand(sy.diff(N,x)*Z-N*sy.diff(Z,x)),x)
    H=H.set_domain(sy.QQ)
    # Count sign-changing roots exactly on the open interval. A rational
    # isolating box touching a boundary does not make its root an endpoint.
    for side,aa,bb in [('heating',sy.Rational(0),sy.Rational(1,2)),('cooling',sy.Rational(1,2),sy.Rational(1))]:
        boundary+=sum(H.eval(t)==0 for t in (aa,bb))
        exact=sign_changes_open(H,aa,bb)>0
        out=classify_mode(np.arange(4)*np.log(2),ww)
        if isinstance(out[side],bool) and out[side]!=exact:errors.append({'w':ww,'side':side,'issue':'exact polynomial disagreement','out':out,'exact':bool(exact)})
        polychecks+=1

# Wronskian identity checked through high-precision differentiation, independent
# of the closed determinant formula. General n verifies the proposed bound basis.
mp.mp.dps=60;wronskian_checks=0
for es in [[0,1,2,3],[0,mp.mpf('.005'),mp.mpf('.916'),6],[0,mp.mpf('.2'),mp.mpf('.31'),1,4]]:
    ee=list(map(mp.mpf,es));nn=len(ee)
    for beta in [mp.mpf('-.8'),mp.mpf('0'),mp.mpf('.7'),mp.mpf('2')]:
        f=[mp.exp(-v*beta) for v in ee];zz=sum(f)
        def pi(i,b):return mp.exp(-ee[i]*b)/sum(mp.exp(-v*b) for v in ee)
        for k in range(1,nn):
            W=mp.det(mp.matrix([[mp.diff(lambda b:pi(j,b),beta,i+1) for j in range(k)] for i in range(k)]))
            vand=mp.mpf(1)
            for i in range(k):
                for j in range(i+1,k):vand*=ee[j]-ee[i]
            tail=sum(f[j]*mp.fprod(ee[j]-ee[i] for i in range(k)) for j in range(k,nn))
            formula=(-1)**(k*(k-1)//2)*vand*mp.fprod(f[:k])*tail/zz**(k+1)
            if abs(W-formula)>mp.mpf('1e-45')*max(abs(formula),mp.mpf('1e-20')):raise ValueError('Wronskian identity failed')
            wronskian_checks+=1
out={'seed':111028,'counts':counts,'exact_polynomial_side_checks':polychecks,'exact_endpoint_zeros_seen':boundary,'wronskian_checks':wronskian_checks,'errors':errors,'examples':examples,'elapsed_seconds':time.time()-start,'scope':'Numerical stress tests, exact-polynomial reference tests, and high-precision identity cross-checks; general proof is separate.'}
(ROOT/'research/classifier_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('FINAL',json.dumps(out),flush=True)
if errors:raise SystemExit('Classifier disagreements require investigation')
