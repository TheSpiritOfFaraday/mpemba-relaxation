"""Independent numerical audit of every published figure CSV row.

Direct exponentiation checks evolution curves; numerical differentiation
checks thermal stationary points. This is numerical correspondence, not
an interval proof of continuum claims. No plotting module is imported.
"""
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/figures'
mp.mp.dps=70
stats={}


def rows(name):
    with (OUT/name).open() as f:
        return list(csv.DictReader(f))


def check(group,actual,expected,tol='2e-11'):
    a,b=mp.mpf(actual),mp.mpf(expected)
    error=abs(a-b)
    if not mp.isfinite(error) or error>mp.mpf(tol):
        raise ArithmeticError((group,str(a),str(b),str(error)))
    record=stats.setdefault(group,{'comparisons':0,'max_absolute_error':mp.mpf(0)})
    record['comparisons']+=1
    record['max_absolute_error']=max(error,record['max_absolute_error'])


def rat(text):
    parts=text.split('/')
    return mp.mpf(parts[0])/(mp.mpf(parts[1]) if len(parts)>1 else 1)


def rational_model(name):
    d=json.loads((ROOT/'research'/name).read_text())
    q=mp.matrix([rat(v) for v in d['q']])
    M=mp.matrix([[rat(v) for v in row] for row in d['generator']])
    pi=q/sum(q)
    initial=[]
    for beta in [mp.mpf(1)/1000,mp.mpf(1)/100]:
        v=mp.matrix([v**beta for v in q]);initial.append(v/sum(v))
    return q,M,pi,initial


def distances(p,pi,weights):
    return (sum(abs(p[i]-pi[i]) for i in range(4))/2,
            sum(p[i]*mp.log(p[i]/pi[i]) for i in range(4)),
            sum((p[i]-pi[i])**2/weights[i] for i in range(4)))


q,M,pi,initial=rational_model('adversarial_certificate.json')
for row in rows('relaxation_data.csv'):
    K=mp.expm(M*mp.mpf(row['time']))
    for j,label in enumerate(['hot','cold']):
        tv,kl,_=distances(K*initial[j],pi,q)
        check('primary direct evolution',row['TV_'+label],tv)
        check('primary direct evolution',row['KL_'+label],kl)
S=mp.matrix([[M[i,j]*mp.sqrt(q[j]/q[i]) for j in range(4)] for i in range(4)])
ev,U=mp.eigsy(S)
w=[U[i,2]/mp.sqrt(q[i]) for i in range(4)];w=[v/w[0] for v in w]
def primary_overlap(beta):
    weights=[v**beta for v in q]
    return mp.fdot(weights,w)/sum(weights)
base=abs(primary_overlap(mp.mpf('.001')))
for row in rows('overlap_data.csv'):
    value=abs(primary_overlap(mp.mpf(row['beta'])))
    check('primary overlap',row['absolute_overlap'],value)
    check('primary overlap',row['enhancement_percent'],100*(value/base-1))

q,M,pi,initial=rational_model('adversarial_rational_mode_certificate.json')
for row in rows('focused_relaxation.csv'):
    K=mp.expm(M*mp.mpf(row['time']))
    values=[distances(K*p,pi,q)[2] for p in initial]
    check('rational direct evolution',row['distance_hot'],values[0])
    check('rational direct evolution',row['distance_cold'],values[1])
    check('rational direct evolution',row['relative_percent'],100*(values[0]/values[1]-1))

# Differentiate the defining thermal sum, independently of the plotted
# pair-expansion derivative and its bisection implementation.
q=[mp.mpf(1),mp.mpf(199)/200,mp.mpf(2)/5,mp.mpf(1)/400]
def mechanism(z,beta):
    w=[1,-mp.mpf(759)/1000,-mp.mpf(48959)/80000-z/160,z]
    weights=[v**beta for v in q]
    return mp.fdot(w,weights)/sum(weights)
zcrit=mp.findroot(lambda z:mp.diff(lambda b:mechanism(z,b),0),(-mp.mpf('.03'),-mp.mpf('.004')))
for row in rows('mechanism_curves.csv'):
    for key,z in [('z_minus025',mp.mpf('-.025')),('critical',zcrit),('z_minus012',mp.mpf('-.012')),('z_minus004',mp.mpf('-.004'))]:
        expected=10000*(abs(mechanism(z,mp.mpf(row['beta'])))-abs(mechanism(z,0)))
        check('mechanism curves',row[key],expected)
for row in rows('mechanism_onset.csv'):
    z,beta=mp.mpf(row['z']),mp.mpf(row['peak_beta'])
    check('mechanism peak derivative',mp.diff(lambda b:mechanism(z,b),beta),0)

def quarter(ell,beta):
    eps=mp.exp(-ell)
    E=[0,-mp.log1p(-eps/3),mp.log(mp.mpf(5)/3),2*ell]
    w=[1,-mp.mpf(4)/5,-mp.mpf(1)/3-4*eps/9+5*eps**3/3,-eps]
    weights=[mp.exp(-beta*e) for e in E]
    return mp.fdot(weights,w)/sum(weights)
for row in rows('quarter_overlap_data.csv'):
    ell=mp.mpf(row['exponent'])*mp.log(10);beta=mp.mpf(row['beta']);hot=1/ell**2
    gain=100*(1-abs(quarter(ell,hot)/quarter(ell,beta)))
    check('quarter thermal curves',row['gain'],gain)
    check('quarter hot preparation',row['hot_beta'],hot)
for index,row in enumerate(rows('quarter_gain_data.csv')):
    ell=(3+mp.mpf(47)*index/100)*mp.log(10);beta=mp.mpf(row['peak_beta'])
    check('quarter maxima',row['energy_span'],2*ell)
    check('quarter maxima',mp.diff(lambda b:quarter(ell,b),beta),0)
    check('quarter maxima',row['maximum_gain'],100*(1-abs(quarter(ell,1/ell**2)/quarter(ell,beta))))

candidate=json.loads((ROOT/'research/physical_rate_refined.json').read_text(),parse_float=Decimal)['result']['parameters']
audit=json.loads((ROOT/'research/physical_rate_refined_audit.json').read_text(),parse_float=Decimal)
if candidate!=audit['parameters']:
    raise ArithmeticError('Physical figure input differs from certified decimal input')
x=list(map(lambda v:mp.mpf(str(v)),candidate));E=[0,x[0],x[0]+x[1],sum(x[:3])];M=mp.matrix(4)
for i in range(4):
    for j in range(i+1,4):
        h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        saddle=E[j]+h
        M[i,j]=mp.exp(-(saddle-E[j]));M[j,i]=mp.exp(-(saddle-E[i]))
for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
def thermal(T):
    v=mp.matrix([mp.exp(-e/T) for e in E]);return v/sum(v)
pi=thermal(1);initial=[thermal(T) for T in x[-2:]]
for row in rows('physical_finite_time.csv'):
    t=mp.mpf(row['time_attempt_units']);K=mp.expm(M*t)
    check('physical time units',t,mp.mpf(row['time_tau'])*mp.mpf(str(audit['nominal']['tau'])))
    for j,label in enumerate(['cold','hot']):
        tv,kl,chi=distances(K*initial[j],pi,pi)
        for key,val in [('KL',kl),('TV',tv),('chi2',chi)]:check('physical direct evolution',row[key+'_'+label],val)
for i,row in enumerate(rows('physical_threshold_savings.csv')):
    tt=audit['nominal']['hitting_times_tau'];rob=audit['robustness']
    check('physical envelope correspondence',row['threshold_nats'],str(rob['thresholds'][i]))
    check('physical envelope correspondence',row['nominal_percent'],str(100*(1-tt[1][i]/tt[0][i])))
    for key,field in [('sample_min_percent','time_saving_envelope_low'),('sample_max_percent','time_saving_envelope_high')]:
        check('physical envelope correspondence',row[key],str(100*rob[field][i]))
receipt={'status':'passed','precision_digits':70,'scope':'Every row of nine figure CSVs checked numerically. Exact claims rely on separate reviewed certificates; sampled envelope uses its separately replayed audit.',
         'groups':{k:{'comparisons':v['comparisons'],'max_absolute_error':str(v['max_absolute_error'])} for k,v in stats.items()},
         'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*sorted(OUT.glob('*.csv'))]}}
(ROOT/'research/detailed_audit/figure_data_checks.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS:',sum(v['comparisons'] for v in stats.values()),'numerical figure-data comparisons')
for key,val in receipt['groups'].items():print(key,val)
