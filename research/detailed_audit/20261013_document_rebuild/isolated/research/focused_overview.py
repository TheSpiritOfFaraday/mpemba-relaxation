"""Reproducible numerical overview; curves illustrate separately proved claims."""
from pathlib import Path
import csv
import json
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/figures'
mp.mp.dps = 70
q = [mp.mpf(1), mp.mpf(199)/200, mp.mpf(2)/5, mp.mpf(1)/400]
w = [mp.mpf(1), -mp.mpf(759)/1000, -mp.mpf(48957)/80000, -mp.mpf(1)/250]
rates = [mp.mpf(1), mp.mpf(28), mp.mpf(249)/10,
         mp.mpf(330154675693)/1470375000, mp.mpf(4983)/200,
         mp.mpf(1527327066)/151990625]
M = mp.matrix(4)
k = 0
for i in range(4):
    for j in range(i+1,4):
        M[i,j] = rates[k]
        M[j,i] = rates[k]*q[j]/q[i]
        k += 1
for j in range(4):
    M[j,j] = -sum(M[i,j] for i in range(4) if i != j)
S = mp.matrix([[M[i,j]*mp.sqrt(q[j]/q[i]) for j in range(4)] for i in range(4)])
ev, U = mp.eigsy(S)
assert abs(ev[2]+mp.mpf(4966671)/250000) < mp.mpf('1e-60')
def gibbs(beta):
    a = [x**beta for x in q]
    return mp.matrix([x/sum(a) for x in a])
bath = mp.matrix([x/sum(q) for x in q])
ph, pc = gibbs(mp.mpf(1)/1000), gibbs(mp.mpf(1)/100)
def spectral(p,t):
    x = mp.matrix([(p[i]-bath[i])/mp.sqrt(q[i]) for i in range(4)])
    coeff = U.T*x
    return mp.matrix([bath[i]+mp.sqrt(q[i])*sum(U[i,k]*mp.exp(ev[k]*t)*coeff[k]
                       for k in range(3)) for i in range(4)])
def distance(p):
    return sum((p[i]-bath[i])**2/q[i] for i in range(4))
max_error = mp.mpf(0)
for t in map(mp.mpf,['0','.03','.1','.2','.35']):
    for p in [ph,pc]:
        a, b = spectral(p,t), mp.expm(M*t)*p
        max_error = max(max_error,max(abs(a[i]-b[i]) for i in range(4)))
assert max_error < mp.mpf('1e-60')
curve=[]
for j in range(181):
    t = mp.mpf(j)/600
    dh,dc=distance(spectral(ph,t)),distance(spectral(pc,t))
    curve.append({'time':float(t),'distance_hot':float(dh),'distance_cold':float(dc),
                  'relative_percent':float(100*(dh/dc-1))})
with (OUT/'focused_relaxation.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=curve[0]);wr.writeheader();wr.writerows(curve)
with (OUT/'mechanism_onset.csv').open() as f:
    onset=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
with (OUT/'quarter_gain_data.csv').open() as f:
    sharp=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titlesize':8.2,'pdf.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(7.0,2.7),layout='constrained')
a=axes[0]
a.plot([r['time'] for r in curve],[r['relative_percent'] for r in curve],color='#176b87',lw=2)
a.axhline(0,color='#444444',lw=.8)
a.fill_between([r['time'] for r in curve],[r['relative_percent'] for r in curve],0,
               where=[r['relative_percent']<0 for r in curve],color='#176b87',alpha=.18)
a.set(title='A  Distance reversal',xlabel=r'Time ($a_{12}^{-1}$)',
      ylabel=r'$100(D_h/D_c-1)$ (%)',xlim=(0,.3))
a.text(.97,.94,'Rational witness',ha='right',va='top',transform=a.transAxes,fontsize=7)
a.text(.97,.22,'Below zero:\nhotter is closer',ha='right',transform=a.transAxes,fontsize=7,color='#176b87')
a=axes[1]
zstar=onset[0]['z']
a.plot([r['z'] for r in onset],[r['peak_beta'] for r in onset],color='#b15e20',lw=2)
a.fill_between([r['z'] for r in onset],0,[r['peak_beta'] for r in onset],color='#b15e20',alpha=.18)
a.axvline(zstar,color='#555555',ls='--',lw=1)
a.set(title='B  Controlled onset',xlabel=r'Mode component $z$',
      ylabel=r'Peak location $\beta_m$',xlim=(-.03,-.004),ylim=(0,.03))
a.set_xticks([-.03,-.02,-.01])
a.text(-.029,.022,'No cooling\nreversal',fontsize=7)
a.text(-.009,.004,'Cooling\nwindow',ha='center',fontsize=7,color='#b15e20')
a=axes[2]
a.plot([r['energy_span'] for r in sharp],[r['maximum_gain'] for r in sharp],color='#7551a6',lw=2)
a.axhline(25,color='#333333',ls='--',lw=1)
a.set(title='C  Sharp amplitude limit',xlabel=r'Energy span $\beta_b(E_4-E_1)$',
      ylabel=r'Amplitude advantage (%)',xlim=(0,240),ylim=(0,28))
a.text(235,25.7,'25% supremum',ha='right',fontsize=7)
a.text(.96,.14,'Finite thermal\npreparations',ha='right',transform=a.transAxes,fontsize=7,color='#7551a6')
fig.savefig(OUT/'focused_overview.pdf')
fig.savefig(OUT/'focused_overview.png',dpi=200)
bh,bc=mp.fdot(w,ph),mp.fdot(w,pc)
receipt={'precision_digits':70,'generator':'Appendix B rational companion',
 'spectral_vs_expm_max_absolute_error':str(max_error),
 'initial_D_hot':str(distance(ph)),'initial_D_cold':str(distance(pc)),
 'amplitude_advantage_percent':str(100*(1-abs(bh)/abs(bc))),
 'asymptotic_distance_advantage_percent':str(100*(1-(bh/bc)**2)),
 'panel_B_source':'output/figures/mechanism_onset.csv',
 'panel_C_source':'output/figures/quarter_gain_data.csv',
 'panel_C_hot_beta':'1/log(1/epsilon)^2; cold at numerical overlap maximum',
 'scope':'Numerical illustration only. No rigorous finite crossing time is inferred.'}
(OUT/'focused_overview_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
