"""Reproducible physical illustrations of the primary exact counterexample.
High-precision mpmath (80 digits); figures are numerical, distinct from exact
certificates. Matrix exponential and symmetric spectral propagation are compared.
"""
from pathlib import Path
import json,csv
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
mp.mp.dps=80
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/figures';OUT.mkdir(exist_ok=True,parents=True)
r=json.loads((ROOT/'research/adversarial_certificate.json').read_text())
def rat(z):
 x=z.split('/');return mp.mpf(x[0])/(mp.mpf(x[1]) if len(x)>1 else 1)
q=mp.matrix([rat(z) for z in r['q']]);M=mp.matrix([[rat(z) for z in row] for row in r['generator']]);pi=q/mp.fsum(q);E=mp.matrix([-mp.log(z) for z in q])
D=mp.diag([mp.sqrt(z) for z in pi]);S=D**-1*M*D
lam,U=mp.eigsy(S);slow=2
w=mp.matrix([U[i,slow]/mp.sqrt(pi[i]) for i in range(4)]);w/=w[0]
def gibbs(beta):
 z=mp.matrix([qi**beta for qi in q]);return z/mp.fsum(z)
def overlap(beta):return mp.fdot(w,gibbs(beta))
bh=mp.mpf(1)/1000;bc=mp.mpf(1)/100
starts=[gibbs(bh),gibbs(bc)];coeff=[U.T*(D**-1*(p-pi)) for p in starts]
def state(t,j):
 # Exclude stationary mode to avoid tiny eigensolver residuals at long times.
 return pi+D*U*mp.matrix([coeff[j][k]*mp.exp(lam[k]*t) if k<3 else 0 for k in range(4)])
def metrics(t,j):
 p=state(mp.mpf(t),j);delta=p-pi
 tv=mp.fsum(abs(x) for x in delta)/2
 kl=mp.fsum(p[i]*mp.log1p(delta[i]/pi[i]) for i in range(4))
 return tv,kl
checks=[]
for t in ['0','.05','.1','.2','.4']:
 for j in range(2):
  direct=mp.expm(M*mp.mpf(t))*starts[j]
  err=max(abs(direct[i]-state(mp.mpf(t),j)[i]) for i in range(4));assert err<mp.mpf('1e-70');checks.append(str(err))
def diff(t,k):return metrics(t,0)[k]-metrics(t,1)[k]
roots=[]
for k in [0,1]:
 scan=np.linspace(0,.8,161);cross=[]
 for x,y in zip(scan[:-1],scan[1:]):
  if diff(str(x),k)*diff(str(y),k)<0:
   l,h=mp.mpf(str(x)),mp.mpf(str(y))
   for _ in range(180):
    mid=(l+h)/2
    if diff(l,k)*diff(mid,k)<=0:h=mid
    else:l=mid
   cross.append((l+h)/2)
 roots.append(cross)
print('crossings TV KL',[[str(z) for z in a] for a in roots],flush=True)
# Scalar b minimum near high-temperature endpoint.
beta_turn=mp.findroot(lambda b:mp.diff(overlap,b),(.01,.05))
rel=1-abs(overlap(bh)/overlap(bc));klrel=1-(overlap(bh)/overlap(bc))**2
shift=mp.log(abs(overlap(bc)/overlap(bh)))/(-lam[slow])
summary={'precision_digits':80,'initial_beta':['1/1000','1/100'],'normalization':'w1=1','lambda2':str(lam[slow]),'overlaps':[str(overlap(bh)),str(overlap(bc))],'relative_amplitude_reduction':str(rel),'asymptotic_KL_reduction':str(klrel),'asymptotic_time_advance':str(shift),'slow_time_constant':str(-1/lam[slow]),'beta_turn_numerical':str(beta_turn),'crossings_numerical':{'TV':[str(z) for z in roots[0]],'KL':[str(z) for z in roots[1]]},'expm_vs_spectral_max_error':str(max(map(mp.mpf,checks))),'samples':{}}
for t in ['0','.05','.1','.15','.2','.25','.3','.4']:
 vals=[metrics(t,j) for j in range(2)];summary['samples'][t]={'hot_TV':str(vals[0][0]),'cold_TV':str(vals[1][0]),'hot_KL':str(vals[0][1]),'cold_KL':str(vals[1][1])}
(OUT/'physics_metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','font.family':'DejaVu Sans'})
colors=['#b33b32','#246ba0']
# Full overlap and enlargement with zero-based enhancement to show small anomaly.
fig,axs=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
bs=np.linspace(0,1,350);ovs=np.array([float(abs(overlap(mp.mpf(str(b))))) for b in bs]);base=float(abs(overlap(bh)))
axs[0].plot(bs,ovs,color='#202833',lw=2);axs[0].set(xlabel=r'Initial inverse temperature $\beta$',ylabel=r'Slow-mode amplitude $|b_2(\beta)|$',title='A  Full cooling interval',xlim=(0,1));axs[0].axvline(1,color='#888',ls=':',lw=1)
zoom=np.linspace(0,.07,250);vals=100*(np.array([float(abs(overlap(mp.mpf(str(b))))) for b in zoom])/base-1)
axs[1].plot(zoom,vals,color='#202833',lw=2);axs[1].axhline(0,color='#999',lw=.8);axs[1].set(xlabel=r'Initial inverse temperature $\beta$',ylabel=r'Amplitude increase from hot state (%)',title='B  High-temperature enlargement',xlim=(0,.07))
for b,c,label in zip([float(bh),float(bc)],colors,[r'$T_h=1000$',r'$T_c=100$']):
 value=float(abs(overlap(mp.mpf(str(b)))));axs[0].scatter(b,value,color=c,s=35,zorder=5);axs[1].scatter(b,100*(value/base-1),color=c,s=40,zorder=5,label=label)
with (OUT/'overlap_data.csv').open('w') as f:
 wr=csv.writer(f);wr.writerow(['beta','absolute_overlap','enhancement_percent'])
 for b,v in zip(bs,ovs):wr.writerow([b,v,100*(v/base-1)])
 for b,v in zip(zoom,vals):wr.writerow([b,base*(1+v/100),v])
axs[1].legend(frameon=False,loc='lower left');axs[1].axvline(float(beta_turn),color='#999',lw=.8,ls=':')
fig.savefig(OUT/'thermal_overlap.pdf');fig.savefig(OUT/'thermal_overlap.png',dpi=220);plt.close(fig)
# Distance values and ratio are both important: main curves almost coincide.
ts=np.linspace(0,.5,350);mt=np.array([[[float(x) for x in metrics(str(t),j)] for j in range(2)] for t in ts])
with (OUT/'relaxation_data.csv').open('w') as f:
 wr=csv.writer(f);wr.writerow(['time','TV_hot','TV_cold','KL_hot','KL_cold'])
 for t,row in zip(ts,mt):wr.writerow([t,row[0,0],row[1,0],row[0,1],row[1,1]])
fig,axs=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
for j in range(2):axs[0].semilogy(ts,mt[:,j,1],color=colors[j],ls='-' if j==0 else '--',label=['Hot initial state','Cool initial state'][j],lw=2)
axs[0].set(xlabel=r'Time $t$ (units $a_{12}^{-1}$)',ylabel=r'$D_{\rm KL}(p(t)\Vert\pi_b)$',title='A  Actual relaxation distance');axs[0].legend(frameon=False)
ratio=100*(mt[:,0,1]/mt[:,1,1]-1);axs[1].plot(ts,ratio,color='#202833',lw=2);axs[1].axhline(0,color='#888',lw=1)
axs[1].axhline(-100*float(klrel),color='#999',lw=1,ls=':',label='Asymptotic limit')
rt=float(roots[1][0]);axs[1].scatter(rt,0,color='#7952a3',zorder=5);axs[1].annotate(f'Crossing at t = {rt:.5f}',(rt,0),xytext=(.23,1.5),arrowprops={'arrowstyle':'->','color':'#555'},fontsize=9)
axs[1].set(xlabel=r'Time $t$ (units $a_{12}^{-1}$)',ylabel=r'$100(D_h/D_c-1)$ (%)',title='B  Sign change resolves the crossing');axs[1].legend(frameon=False,loc='upper right')
fig.savefig(OUT/'relaxation_crossing.pdf');fig.savefig(OUT/'relaxation_crossing.png',dpi=220);plt.close(fig)
print(json.dumps({k:v for k,v in summary.items() if k not in ['samples']},indent=2),flush=True)
