"""Numerical illustrations for the analytically certified quarter-limit family.

These curves are explanatory, not acceptance certificates. All overlaps and
stationary points are evaluated with 80-digit arithmetic. Initial hot beta is
L^-2, L=log(1/epsilon); the curve therefore uses finite positive preparations.
"""
import csv,json
from pathlib import Path
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp.mp.dps=80
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/figures';OUT.mkdir(exist_ok=True)

def system(L):
    e=mp.exp(-L)
    q=[mp.mpf(1),1-e/3,mp.mpf(3)/5,e*e]
    w=[mp.mpf(1),-mp.mpf(4)/5,-mp.mpf(1)/3-4*e/9+5*e**3/3,-e]
    E=[-mp.log(v) for v in q]
    if abs(sum(a*b for a,b in zip(q,w)))>mp.mpf('1e-70'):
        raise ArithmeticError('bath orthogonality')
    def values(beta):
        f=[mp.exp(-beta*v) for v in E];Z=sum(f)
        b=sum(fi*wi for fi,wi in zip(f,w))/Z
        mean=sum(fi*Ei for fi,Ei in zip(f,E))/Z
        der=-sum(fi*wi*(Ei-mean) for fi,wi,Ei in zip(f,w,E))/Z
        return b,der
    return values

def peak(L):
    values=system(L);hot=1/L**2;bh=values(hot)[0]
    if values(hot)[1]>=0:return hot,mp.mpf(0)
    lo=hot;hi=mp.mpf(1)
    for _ in range(240):
        mid=(lo+hi)/2
        if values(mid)[1]<0:lo=mid
        else:hi=mid
    beta=(lo+hi)/2;bc=values(beta)[0]
    return beta,100*(1-abs(bh)/abs(bc))

rows=[]
for exponent in [3,10,50]:
    L=exponent*mp.log(10);values=system(L);hot=1/L**2;bh=values(hot)[0]
    beta_peak,gain=peak(L)
    for j in range(151):
        beta=hot+(mp.mpf('0.35')-hot)*j/150
        b=values(beta)[0]
        rows.append(dict(exponent=exponent,beta=float(beta),gain=float(100*(1-abs(bh)/abs(b))),
            hot_beta=float(hot),peak_beta=float(beta_peak),peak_gain=float(gain)))
with (OUT/'quarter_overlap_data.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)

maxima=[]
for j in range(101):
    L=(3+mp.mpf(47)*j/100)*mp.log(10)
    beta,gain=peak(L)
    maxima.append(dict(energy_span=float(2*L),maximum_gain=float(gain),peak_beta=float(beta)))
with (OUT/'quarter_gain_data.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(maxima[0]));wr.writeheader();wr.writerows(maxima)

fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for exponent in [3,10,50]:
    rr=[r for r in rows if r['exponent']==exponent]
    axes[0].plot([r['beta'] for r in rr],[r['gain'] for r in rr],label=rf'$\epsilon=10^{{-{exponent}}}$')
axes[0].set(xlim=(0,.35),ylim=(0,26),xlabel=r'Initial inverse temperature $\beta$',ylabel='Amplitude advantage (%)')
axes[0].legend(frameon=False)
axes[1].plot([r['energy_span'] for r in maxima],[r['maximum_gain'] for r in maxima],color='black')
axes[1].axhline(25,color='gray',linestyle='--',label='Proved limiting advantage')
axes[1].set(xlabel=r'Bath-scaled energy span $\beta_b(E_4-E_1)$',ylabel='Maximum amplitude advantage (%)',ylim=(0,26))
axes[1].legend(frameon=False,loc='lower right')
fig.savefig(OUT/'quarter_limit.png',dpi=180)
fig.savefig(OUT/'quarter_limit.pdf')
receipt=dict(precision_digits=80,thermal_points=len(rows),maxima_points=len(maxima),
    hot_preparation='beta_h=1/log(1/epsilon)^2',
    representative_peaks=[{k:r[k] for k in ['exponent','hot_beta','peak_beta','peak_gain']} for r in rows[::151]],
    scope='Numerical illustration of the unperturbed certified family, not additional rigorous acceptance')
(OUT/'quarter_plot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
