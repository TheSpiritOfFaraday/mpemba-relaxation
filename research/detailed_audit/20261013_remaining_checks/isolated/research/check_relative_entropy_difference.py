"""Independent direct-matrix audit of the cancellation-preserving KL identity."""
import itertools
import json
import subprocess
from pathlib import Path
import mpmath as mp
import numpy as np
from check_physical_box_spectral import direct,enclosed
ROOT=Path(__file__).resolve().parent
mp.mp.dps=70

def difference(x,u):
    _,D=direct(x,u);_,H=direct(x,mp.mpf('.95')*u)
    return D[0]-H[1]

if __name__=='__main__':
    x=[mp.mpf(str(v)) for v in json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']]
    rng=np.random.default_rng(12102044)
    signs=list(itertools.product([-1,1],repeat=8))+list(rng.uniform(-1,1,(24,8)))
    points=[x]+[[int(k>=6)+(a-int(k>=6))*(1+mp.mpf('.05')*mp.mpf(str(s[k]))) for k,a in enumerate(x)] for s in signs]
    checks=gradients=0;maxwidth=0.;step=mp.mpf('1e-20')
    for n,y in enumerate(points):
        for ustr in ['.1','.2','.4','1.5']:
            u=mp.mpf(ustr)
            r=json.loads(subprocess.check_output([str(ROOT.parent/'tmp/check_relative_entropy_difference'),ustr]+list(map(str,y)),text=True))
            F=difference(y,u)
            assert enclosed(F,r['value']),(n,ustr,str(F),r['value'])
            checks+=1;maxwidth=max(maxwidth,r['value'][1]-r['value'][0])
            if n in [0,1,128,256,280] and ustr in ['.2','1.5']:
                for k in range(9):
                    yp=y.copy();ym=y.copy();up=u;um=u
                    if k<8:yp[k]+=step;ym[k]-=step
                    else:up+=step;um-=step
                    d=(difference(yp,up)-difference(ym,um))/(2*step)
                    assert enclosed(d,r['gradient'][k]),(n,ustr,k,str(d),r['gradient'][k])
                    gradients+=1
        if n%50==0:print('audited points',n+1,flush=True)
    receipt=dict(precision_digits=70,points=len(points),value_enclosure_checks=checks,derivative_enclosure_checks=gradients,
                 maximum_value_interval_width=maxwidth,all_passed=True,seed=12102044,
                 method='Independent direct 4x4 matrix exponentials; derivative central differences with step 1e-20',
                 scope='Point audit of the KL difference identity and its derivatives; not a full-box certificate')
    (ROOT/'relative_entropy_difference_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
