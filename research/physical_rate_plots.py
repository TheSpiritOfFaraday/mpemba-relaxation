"""Plot the refined candidate; sampled envelopes are explicitly not certificates."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from physical_rate_model import Evolution

ROOT=Path(__file__).resolve().parents[1]
audit=json.loads((ROOT/'research/physical_rate_refined_audit.json').read_text())
e=Evolution(audit['parameters'])
s=np.linspace(0,1.5,151)
D=e.distances(s)
tv=e.distances(s,'tv')
chi=e.distances(s,'chi')
out=ROOT/'output/figures'
with (out/'physical_finite_time.csv').open('w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['time_tau','time_attempt_units','KL_cold','KL_hot','TV_cold','TV_hot','chi2_cold','chi2_hot'])
    for j,t in enumerate(s):
        w.writerow([t,t*e.tau,*D[:,j],*tv[:,j],*chi[:,j]])
rob=audit['robustness']
delta=np.array(rob['thresholds'])
tt=np.array(audit['nominal']['hitting_times_tau'])
gain=100*(1-tt[1]/tt[0])
low=100*np.array(rob['time_saving_envelope_low'])
high=100*np.array(rob['time_saving_envelope_high'])
with (out/'physical_threshold_savings.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['threshold_nats','nominal_percent','sample_min_percent','sample_max_percent'])
    w.writerows(zip(delta,gain,low,high))
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(8,3.2),layout='constrained')
ax[0].semilogy(s,D[0],label='Cold preparation',color='#276b8e',lw=2)
ax[0].semilogy(s,D[1],label='Hot preparation',color='#c26026',lw=2)
ax[0].axvspan(.4,.8,color='#777777',alpha=.12,label='Certified distance interval')
ax[0].axhline(.001,color='#777777',ls=':',lw=1)
ax[0].set(xlabel=r'Time / $\tau$',ylabel='KL distance (nats)',title='A  Finite-time relaxation',xlim=(0,1.5),ylim=(.0006,1))
ax[0].legend(frameon=False,fontsize=7,loc='upper right')
ax[1].fill_between(delta,low,high,color='#a7c2d0',alpha=.55,label='Sampled ±5% envelope')
ax[1].semilogx(delta,gain,color='#276b8e',lw=2,label='Nominal')
ax[1].axhline(5,color='#a14922',ls='--',lw=1,label='Certified 5% floor (full box)')
ax[1].set(xlabel='Common KL threshold (nats)',ylabel='Time saved (%)',title='B  Cooling-time savings',xlim=(.001,.01),ylim=(0,max(high)+3))
ax[1].legend(frameon=False,fontsize=7,loc='upper left')
ax[1].set_xticks([.001,.003,.01],labels=[r'$10^{-3}$',r'$3\times10^{-3}$',r'$10^{-2}$'])
ax[1].xaxis.set_minor_formatter(NullFormatter())
fig.savefig(out/'physical_finite_time.pdf')
fig.savefig(out/'physical_finite_time.png',dpi=200)

def coords(x,y):
    return ' '.join(f'({a:.12g},{b:.12g})' for a,b in zip(x,y))

tex=r'''\begin{figure}[htbp]
\centering
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=1.3cm},
width=0.45\textwidth,height=6.2cm,tick label style={font=\scriptsize},
label style={font=\small},title style={font=\small},
legend style={font=\scriptsize,draw=none}]
\nextgroupplot[title={A: Finite-time relaxation},xlabel={$t/\tau$},
ylabel={KL distance (nats)},ymode=log,xmin=0,xmax=1.5,ymin=0.0006,ymax=1,
legend pos=north east]
\addplot[fill=gray!15,draw=none,forget plot] coordinates {(.4,.0006) (.8,.0006) (.8,1) (.4,1)} \closedcycle;
'''
tex+=r'\addplot[blue,thick] coordinates {'+coords(s,D[0])+r'};\addlegendentry{Cold}'+'\n'
tex+=r'\addplot[orange,thick] coordinates {'+coords(s,D[1])+r'};\addlegendentry{Hot}'+'\n'
tex+=r'''\nextgroupplot[title={B: Cooling-time savings},xlabel={KL threshold (nats)},
ylabel={Time saved (\%)},xmode=log,xmin=.001,xmax=.01,ymin=0,ymax=60,
legend pos=north west]
'''
tex+=r'\addplot[fill=blue!15,draw=none,area legend] coordinates {'+coords(delta,low)+' '+coords(delta[::-1],high[::-1])+r'} \closedcycle;\addlegendentry{Sampled $\pm5\%$ range}'+'\n'
tex+=r'\addplot[blue,thick] coordinates {'+coords(delta,gain)+r'};\addlegendentry{Nominal}'+'\n'
tex+=r'''\addplot[orange,dashed] coordinates {(.001,5) (.01,5)};\addlegendentry{Certified $5\%$ floor}
\end{groupplot}
\end{tikzpicture}
\caption{The constrained Arrhenius candidate of Section~\ref{sec:physicalresults}.
The shaded time interval in A has an interval certificate of at least
$10\%$ KL-distance advantage with cold distance above $10^{-3}$ nats.
The distance claim now holds throughout the simultaneous $\pm5\%$ box.
In B, the nominal model has a certified lower bound of $8.12\%$ time savings
throughout the threshold band. The shaded parameter envelope uses 256 corners
and 3,000 sampled interior points at 31 thresholds; it is not a certificate
for the whole parameter box. The dashed $5\%$ floor is separately
certified throughout the full simultaneous $\pm5\%$ box for every
threshold in the band (Section~\ref{sec:physicalboxcertificate}).
The approximate nominal relaxation time is
$\tau=147.232\,\nu_0^{-1}$.}
\label{fig:physicalfinite}
\end{figure}
'''
(ROOT/'research/physical_finite_time_figure.tex').write_text(tex)
print('Wrote physical finite-time figure, inline TeX and inspectable CSV data.')
