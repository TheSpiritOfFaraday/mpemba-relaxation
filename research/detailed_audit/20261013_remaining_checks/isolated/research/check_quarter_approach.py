"""Check the approach-rate algebra and independently audit quarter-plot data.

Symbolic identities support the derivation; finite high-precision residuals
illustrate, but do not prove, asymptotic big-O bounds or uniqueness.
"""
import csv
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 100
checks = []


def require(name, condition):
    if not condition:
        raise ArithmeticError(name)
    checks.append(name)


a, L, y, x = s.symbols('a L y x', positive=True)
R = (y / 3 - s.Rational(1, 5)) / (2 + y + x)
derivative = -a*y*s.diff(R, y) - 2*L*x*s.diff(R, x)
numerator = x*((10*L-5*a)*y-6*L)-13*a*y
require('stationarity numerator identity', s.cancel(15*(2+y+x)**2*derivative-numerator) == 0)
require('reference amplitude', R.subs({y: 1, x: 0}) == s.Rational(2, 45))
require('relative beta Taylor coefficient', s.simplify((-a*y*s.diff(R,y)/R).subs({y:1,x:0})) == -13*a/6)
require('relative tail Taylor coefficient', s.simplify((s.diff(R,x)/R).subs({y:1,x:0})) == -s.Rational(1,3))
require('infinite temperature ratio', s.Rational(1,30)/s.Rational(2,45) == s.Rational(3,4))
z=s.symbols('z')
require('gain expansion coefficient', s.diff(1-s.Rational(3,4)/(1-z),z).subs(z,0) == -s.Rational(3,4))


def model(ell):
    eps = mp.exp(-ell)
    energies = [mp.mpf(0), -mp.log1p(-eps/3), mp.log(mp.mpf(5)/3), 2*ell]
    mode = [mp.mpf(1), -mp.mpf(4)/5, -mp.mpf(1)/3-4*eps/9+5*eps**3/3, -eps]

    def overlap(beta):
        weights = [mp.exp(-beta*energy) for energy in energies]
        return mp.fsum(w*f for w,f in zip(weights,mode))/mp.fsum(weights)

    # Numerical differentiation is independent of the covariance formula
    # used by the plotting script.
    def slope(beta):
        return mp.diff(overlap, beta)

    lo, hi = mp.mpf(0), mp.mpf(1)
    require('endpoint signs L='+mp.nstr(ell,8), slope(lo)<0<slope(hi))
    for _ in range(250):
        mid=(lo+hi)/2
        if slope(mid)<0:
            lo=mid
        else:
            hi=mid
    beta=(lo+hi)/2
    return overlap, beta


records=[]
aa=mp.log(mp.mpf(5)/3)
for ell in map(mp.mpf, [10,30,100,300,1000,3000,10000]):
    overlap,beta=model(ell)
    logterm=mp.log(4*ell/(13*aa))
    predicted_beta=logterm/(2*ell)
    gain=1-abs(overlap(0)/overlap(beta))
    predicted_gain=mp.mpf(1)/4-13*aa/(16*ell)*(logterm+1)
    records.append({k:mp.nstr(v,35) for k,v in dict(
        L=ell, peak_beta=beta, gain_infinite_hot_reference=gain,
        beta_leading=predicted_beta, gain_leading=predicted_gain,
        scaled_beta_error=(beta-predicted_beta)*ell**2/mp.log(ell),
        scaled_gain_error=(gain-predicted_gain)*ell**2/mp.log(ell)**2).items()})

plot_errors=[]
with (ROOT/'output/figures/quarter_gain_data.csv').open() as f:
    rows=list(csv.DictReader(f))
for index in [0,25,50,75,100]:
    ell=(3+mp.mpf(47)*index/100)*mp.log(10)
    overlap,beta=model(ell)
    gain=100*(1-abs(overlap(1/ell**2)/overlap(beta)))
    row=rows[index]
    err=max(abs(mp.mpf(row['peak_beta'])-beta), abs(mp.mpf(row['maximum_gain'])-gain), abs(mp.mpf(row['energy_span'])-2*ell))
    require('independent plotted maximum row '+str(index), err<mp.mpf('5e-13'))
    plot_errors.append({'row':index,'max_decimal_difference':mp.nstr(err,20)})

sources=['research/quarter_limit_family.tex','research/quarter_limit_plots.py','research/check_quarter_approach.py','output/figures/quarter_gain_data.csv']
receipt=dict(checks=checks, symbolic_identity_count=6, numerical_precision_digits=100,
             asymptotic_diagnostics=records, plot_comparisons=plot_errors,
             scope='Six exact algebra identities; finite numerical diagnostics and five independently recomputed plotted maxima. Analytic remainder and uniqueness arguments reviewed separately.',
             source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources})
out=ROOT/'research/detailed_audit/quarter_approach_checks.json'
out.write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS:',len(checks),'checks; six symbolic identities; five plotted maxima')
for row in records:
    print(row['L'], row['scaled_beta_error'], row['scaled_gain_error'])
