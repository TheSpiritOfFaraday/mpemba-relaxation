"""Numerical illustration and independent checks of the certified family.
Floating-point plots are illustrative; mechanism_certificate.py proves claims.
"""
from pathlib import Path
import csv,json
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/figures';OUT.mkdir(exist_ok=True)
mp.mp.dps=65
q=list(map(mp.mpf,['1','.995','.4','.0025']));e=[-mp.log(v) for v in q]
p=mp.mpf(759)/1000;r0=mp.mpf(48959)/80000;eta=mp.mpf(1)/160
A=(1-p-r0)*sum(e)+4*(p*e[1]+r0*e[2]);B=(3+eta)*e[3]-(1-eta)*e[1]-(1+3*eta)*e[2];crit=A/B

def mode(z):return [mp.mpf(1),-p,-r0-eta*z,z]
def overlap(z,b):
    v=[x**b for x in q];return mp.fdot(mode(z),v)/sum(v)
def D(z,b):
    w=mode(z)
    return sum((w[i]-w[j])*(e[j]-e[i])*mp.exp(-b*(e[i]+e[j])) for i in range(4) for j in range(i+1,4))
def peak(z):
    if z<=crit:return mp.mpf(0)
    l=mp.mpf(0);h=mp.mpf(1)
    for _ in range(150):
        m=(l+h)/2
        if D(z,m)<0:l=m
        else:h=m
    return (l+h)/2

audit=[]
# Compare exact mode predictions with separately assembled numeric symmetric matrices.
for zz in ['-.03','-.025','-.02','-.018','-.017','-.015','-.01','-.004']:
    z=mp.mpf(zz);wf=mode(z)
    rates=[mp.mpf(1),mp.mpf(28),mp.mpf('24.9'),(528238451563-2257386450*z)/(200000*(11761-500*z)),mp.mpf('24.915'),2*(1550*z*z-5989657*z+1197903)/(5*(80500*z+48959))]
    M=mp.matrix(4)
    k=0
    for i in range(4):
        for j in range(i+1,4):M[i,j]=rates[k];M[j,i]=rates[k]*q[j]/q[i];k+=1
    for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4) if i!=j)
    S=mp.matrix([[M[i,j]*mp.sqrt(q[j]/q[i]) for j in range(4)] for i in range(4)])
    vals,U=mp.eigsy(S);pred=-(1550*z+3973343)/200000
    if abs(vals[2]-pred)>mp.mpf('1e-55'):raise ValueError('Independent spectral check')
    left=[U[i,2]/mp.sqrt(q[i]) for i in range(4)];left=[v/left[0] for v in left]
    if max(abs(left[i]-wf[i]) for i in range(4))>mp.mpf('1e-55'):raise ValueError('Independent mode check')
    # At a nonzero beta use numerical differentiation rather than the pair expansion.
    if abs(mp.diff(lambda b:overlap(z,b),mp.mpf('.03'))-D(z,mp.mpf('.03'))/sum(x**mp.mpf('.03') for x in q)**2)>mp.mpf('1e-55'):raise ValueError('Derivative check')
    bm=peak(z);gain=abs(overlap(z,bm))-abs(overlap(z,0))
    audit.append({'z':zz,'peak_beta':str(bm) if z>crit else None,'gain_over_infinite_temperature_limit':str(gain),'slow_eigenvalue':str(pred)})

curves=[];zs=[mp.mpf('-.025'),crit,mp.mpf('-.012'),mp.mpf('-.004')]
for b in np.linspace(0,.05,151):
    beta=mp.mpf(str(b));curves.append([b]+[float(10000*(abs(overlap(z,beta))-abs(overlap(z,0)))) for z in zs])
phase=[]
for zv in np.linspace(float(crit),-.004,101):phase.append([zv,float(peak(mp.mpf(str(zv))))])
with (OUT/'mechanism_curves.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['beta','z_minus025','critical','z_minus012','z_minus004']);w.writerows(curves)
with (OUT/'mechanism_onset.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['z','peak_beta']);w.writerows(phase)
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
arr=np.array(curves)
for i,label in enumerate(['z = −0.025','z = z*','z = −0.012','z = −0.004']):ax[0].plot(arr[:,0],arr[:,i+1],label=label)
ax[0].axhline(0,color='gray',linewidth=.6);ax[0].set(xlabel='Initial inverse temperature β',ylabel='10⁴ (|b(β)| − |b(0)|)',title='A. Emergence of a thermal maximum');ax[0].legend(fontsize=8)
arr=np.array(phase);ax[1].plot(arr[:,0],arr[:,1],color='black');ax[1].fill_between(arr[:,0],0,arr[:,1],alpha=.18,color='tab:orange',label='|b| increases with β')
ax[1].axvline(float(crit),color='gray',linestyle='--');ax[1].set(xlim=(-.03,-.004),ylim=(0,.027),xlabel='Controlled highest-state mode value z',ylabel='Turning inverse temperature βmax',title='B. Predictive onset');ax[1].legend(fontsize=8)
fig.savefig(OUT/'mechanism_onset.pdf');fig.savefig(OUT/'mechanism_onset.png',dpi=170)
(ROOT/'research/mechanism_numeric_audit.json').write_text(json.dumps({'precision_digits':65,'comparison_tolerance':'1e-55','checks':'8 symmetric eigensystem and differentiated-overlap comparisons passed','onset':str(crit),'sample_points':audit,'scope':'Numerical cross-checks and illustrative plots; not the continuous-parameter proof'},indent=2)+'\n')
print('Eight 65-digit independent spectral/mode/derivative checks passed; figures and CSVs written.')
