"""Explicit short-lived activated intermediates for the four-well network.

For each undirected edge ij add an intermediate s_ij at saddle energy B_ij
with phase-space weight eta. Well weights are one. Bath transitions are
well -> s: 2 exp(-(B_ij-Ewell)); s -> either well: 2/eta.
Detailed balance holds for weights g exp(-E/T_b). Schur elimination gives
the original well generator exactly in the trace clock; laboratory time
retains the intermediate residence time. Full preparations are Gibbs on
all ten states, including intermediates, at each preparation temperature.
This is a mesoscopic kinetic extension, not a derived continuum potential.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from physical_rate_model import construct

ROOT=Path(__file__).resolve().parent


class ActivatedEvolution:
    def __init__(self,parameters,eta):
        E,B,coarse=construct(parameters)
        pairs=list(itertools.combinations(range(4),2))
        self.E=np.r_[E,[B[i,j] for i,j in pairs]]
        self.weights=np.r_[np.ones(4),np.full(6,eta)]
        M=np.zeros((10,10))
        for k,(i,j) in enumerate(pairs,4):
            for a in [i,j]:
                M[k,a]=2*np.exp(-(self.E[k]-self.E[a]))
                M[a,k]=2/eta
        M-=np.diag(M.sum(axis=0))
        self.M=M
        self.pi=self.thermal(1)
        root=np.sqrt(self.pi)
        self.root=root
        S=M*root[None,:]/root[:,None]
        self.lam,self.U=eigh(S)
        self.tau=-1/self.lam[-2]
        self.lam[-1]=0
        self.initial=np.array([self.thermal(T) for T in parameters[-2:]])
        self.coeff=((self.initial-self.pi)/root)@self.U
        self.coeff[:,-1]=0
        self.trace=M[:4,:4]-M[:4,4:]@np.linalg.solve(M[4:,4:],M[4:,:4])
        assert np.max(np.abs(self.trace-coarse))<1e-14
        self.balance_error=np.max(np.abs(M*self.pi[None,:]-M.T*self.pi[:,None]))

    def thermal(self,T):
        w=self.weights*np.exp(-self.E/T)
        return w/sum(w)

    def distances(self,t):
        t=np.atleast_1d(t)
        dev=np.einsum('ik,tk,jk,j->itj',self.coeff,np.exp(t[:,None]*self.lam),self.U,self.root)
        z=dev/self.pi
        return np.sum(self.pi*((1+z)*np.log1p(z)-z),axis=-1)

    def hit(self,delta,j):
        return brentq(lambda t:self.distances(t)[j,0]-delta,0,20*self.tau)

    def metrics(self):
        thresholds=np.geomspace(.001,.01,41)
        times=np.array([[self.hit(d,j) for d in thresholds] for j in range(2)])
        s=np.linspace(.01,5,1000)
        D=self.distances(s*self.tau)
        ok=(D[1]<=.9*D[0])&(D[0]>=.001)
        ids=np.flatnonzero(ok)
        blocks=np.split(ids,np.flatnonzero(np.diff(ids)!=1)+1)
        block=max(blocks,key=len)
        return dict(initial_KL=self.distances(0)[:,0].tolist(),tau=float(self.tau),
                    eigenvalues=self.lam.tolist(),saddle_population_initial=self.initial[:,4:].sum(axis=1).tolist(),
                    saddle_population_bath=float(self.pi[4:].sum()),
                    minimum_time_saving_sampled=float(min(1-times[1]/times[0])),
                    distance_window_tau_sampled=None if not len(block) else [float(s[block[0]]),float(s[block[-1]])],
                    thresholds=thresholds.tolist(),times_attempt_units=times.tolist(),
                    detailed_balance_residual=float(self.balance_error))


if __name__=='__main__':
    parameters=json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']
    results=[]
    for eta in [.01,.025,.05,.1,.2,1.]:
        m=ActivatedEvolution(parameters,eta).metrics()
        results.append(dict(eta=eta,**m))
        print(eta,m['minimum_time_saving_sampled'],m['saddle_population_initial'],flush=True)
    receipt=dict(model="four wells with six explicit activated intermediates",
                 parameters=parameters,clock="physical time including intermediate residence",
                 preparation="full ten-state Gibbs distribution, g_well=1, g_intermediate=eta",
                 physical_scope="phenomenological kinetic extension; not continuum or experiment",
                 results=results,certified=False)
    (ROOT/'activated_intermediate_results.json').write_text(json.dumps(receipt,indent=2)+'\n')
