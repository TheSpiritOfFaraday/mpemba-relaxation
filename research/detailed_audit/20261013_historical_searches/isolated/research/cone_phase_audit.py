"""Classify ordinary cooling/heating in triplet-failing ray mixtures.
Numerical exploration only; no claim of exhaustive rate coverage.
"""
import json,time
from pathlib import Path
from collections import Counter
import numpy as np
from cone_search import rays,generator,pairs,constraints
from four_state_classifier import classify_mode
rng=np.random.default_rng(111029);start=time.time();counts=Counter();examples={};invalid=0
for case in range(1000):
    E=np.r_[0,np.cumsum(10**rng.uniform(-3,.85,3))];q=np.exp(-E);rs=rays(q)
    if len(rs)==0:continue
    aa=np.exp(rng.uniform(-10,4,(500,len(rs))))@rs
    reset=np.array([q[i] for i,j in pairs]);rho=np.maximum(0,(-aa/reset).max(axis=1))+1
    aa+=rho[:,None]*reset
    MM=np.array([generator(q,a) for a in aa]);lam,V=np.linalg.eigh(MM*np.sqrt(q[None,None,:]/q[None,:,None]))
    for j,a in enumerate(aa):
        if (constraints(q)@a).max()>1e-6*max(a):invalid+=1;continue
        if lam[j,-2]-lam[j,-3]<1e-8*max(abs(lam[j])):counts['near_degenerate_skipped']+=1;continue
        w=V[j,:,-2]/np.sqrt(q);w*=1 if w[0]>0 else -1
        result=classify_mode(E,w);key=str((result['cooling'],result['heating']));counts[key]+=1
        if key not in examples or (result['heating'] is True and len(examples)<20):
            examples[key if key not in examples else key+str(case)]=dict(E=E.tolist(),rates=a.tolist(),w=w.tolist(),eigenvalues=lam[j].tolist(),classification=result)
    if case%100==0:print('progress',case,dict(counts),flush=True)
out=dict(seed=111029,counts=dict(counts),invalid=invalid,examples=examples,elapsed_seconds=time.time()-start,scope='Numerical ray-mixture sample; finite simple-mode classifier')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
