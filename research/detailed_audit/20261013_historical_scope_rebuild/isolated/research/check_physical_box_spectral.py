"""Independent 70-digit direct-exponential audit of spectral interval code.

Point/derivative checks supplement, but do not replace, the enclosure proof.
The coverage executable alone determines whether a whole box was covered.
"""
import itertools
import json
import subprocess
from pathlib import Path
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parent
mp.mp.dps=70

def direct(x,u):
    E=[mp.mpf(0),x[0],x[0]+x[1],sum(x[:3])]
    M=mp.matrix(4)
    for i,j in itertools.combinations(range(4),2):
        h=x[3] if (i,j)==(0,3) else x[4] if (i,j)==(1,2) else x[5]
        M[i,j]=mp.exp(-h);M[j,i]=mp.exp(-h-E[j]+E[i])
    for j in range(4):M[j,j]=-sum(M[i,j] for i in range(4))
    def thermal(T):
        q=mp.matrix([mp.exp(-e/T) for e in E]);return q/sum(q)
    pi=thermal(1);S=mp.matrix(4)
    for i in range(4):
        for j in range(4):S[i,j]=M[i,j]*mp.sqrt(pi[j]/pi[i])
    eig=mp.eigsy(S,eigvals_only=True);slow=eig[2]
    P=mp.expm(M*(-u/slow));D=[]
    for T in x[-2:]:
        p=P*thermal(T);D.append(sum(p[i]*mp.log(p[i]/pi[i]) for i in range(4)))
    return slow,D

def enclosed(value,pair):
    return mp.mpf(pair[0])<=value<=mp.mpf(pair[1])

if __name__=='__main__':
    nominal=json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']
    x=[mp.mpf(str(v)) for v in nominal]
    rng=np.random.default_rng(12102042)
    signs=[list(s) for s in itertools.product([-1,1],repeat=8)]+list(rng.uniform(-1,1,(24,8)))
    points=[x]
    for s in signs:
        y=[]
        for k,a in enumerate(x):
            offset=int(k>=6);y.append(offset+(a-offset)*(1+mp.mpf('.05')*mp.mpf(str(s[k]))))
        points.append(y)
    checks=0;gradient_checks=0;maxwidth=0.;h=mp.mpf('1e-20')
    for n,y in enumerate(points):
        for ustr in ['.02','.4','1.5','3']:
            u=mp.mpf(ustr)
            args=[str(ROOT.parent/'tmp/check_physical_box_spectral'),ustr]+[str(a) for a in y]
            result=json.loads(subprocess.check_output(args,text=True))
            slow,D=direct(y,u)
            assert enclosed(slow,result['slow']),(n,ustr,'slow')
            for j in range(2):
                bounds=result['KL'][j]['value']
                assert enclosed(D[j],bounds),(n,ustr,j,str(D[j]),bounds)
                maxwidth=max(maxwidth,bounds[1]-bounds[0]);checks+=1
            if n in [0,1,128,256,280] and ustr in ['.4','1.5']:
                for k in range(9):
                    yp=y.copy();ym=y.copy();up=u;um=u
                    if k<8:yp[k]+=h;ym[k]-=h
                    else:up+=h;um-=h
                    _,dp=direct(yp,up);_,dm=direct(ym,um)
                    for j in range(2):
                        derivative=(dp[j]-dm[j])/(2*h)
                        assert enclosed(derivative,result['KL'][j]['gradient'][k]),(n,ustr,j,k,str(derivative),result['KL'][j]['gradient'][k])
                        gradient_checks+=1
        if n%50==0:print('audited points',n+1,flush=True)
    receipt=dict(precision_digits=70,method='independent direct 4x4 matrix exponential and symmetric eigensolve',
                 parameter_points=len(points),normalized_times=['.02','.4','1.5','3'],
                 distance_enclosure_checks=checks,derivative_enclosure_checks=gradient_checks,
                 derivative_method='central difference with step 1e-20; numerical cross-check, not a derivative certificate',
                 maximum_point_KL_interval_width=maxwidth,all_checks_passed=True,
                 seed=12102042,point_audit_is_not_a_box_certificate=True)
    (ROOT/'physical_box_spectral_independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
