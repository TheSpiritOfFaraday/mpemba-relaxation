"""Independent 70-digit direct-exponential checks of temperature derivatives."""
import itertools
import json
import subprocess
from pathlib import Path
import mpmath as mp
import numpy as np
from check_physical_box_spectral import direct,enclosed
ROOT=Path(__file__).resolve().parent
mp.mp.dps=70
TEMP_STEP=mp.mpf('1e-14')

def partials(x,u,b):
    _,D=direct(x,u);xp=x.copy();xm=x.copy()
    xp[6+b]+=TEMP_STEP;xm[6+b]-=TEMP_STEP
    _,plus=direct(xp,u);_,minus=direct(xm,u)
    return [D[b],(plus[b]-minus[b])/(2*TEMP_STEP),(plus[b]-2*D[b]+minus[b])/TEMP_STEP**2]

if __name__=='__main__':
    x=[mp.mpf(str(v)) for v in json.loads((ROOT/'physical_rate_refined.json').read_text())['result']['parameters']]
    rng=np.random.default_rng(12102045)
    signs=list(itertools.product([-1,1],repeat=8))+list(rng.uniform(-1,1,(24,8)))
    points=[x]+[[int(k>=6)+(a-int(k>=6))*(1+mp.mpf('.05')*mp.mpf(str(s[k]))) for k,a in enumerate(x)] for s in signs]
    checked=gradients=0;step=mp.mpf('1e-12')
    for n,y in enumerate(points):
        for ustr in ['.1','.23','.4','1.5']:
            u=mp.mpf(ustr)
            result=json.loads(subprocess.check_output([str(ROOT.parent/'tmp/check_temperature_reduction'),ustr]+list(map(str,y)),text=True))
            for b in range(2):
                values=partials(y,u,b)
                for order,value in enumerate(values):
                    assert enclosed(value,result[b][order]['value']),(n,ustr,b,order,str(value),result[b][order]['value'])
                    checked+=1
                if n in [0,1,128,256,280] and ustr in ['.23','1.5']:
                    for k in range(9):
                        yp=y.copy();ym=y.copy();up=u;um=u
                        if k<8:yp[k]+=step;ym[k]-=step
                        else:up+=step;um-=step
                        dp=partials(yp,up,b);dm=partials(ym,um,b)
                        for order in [1,2]:
                            derivative=(dp[order]-dm[order])/(2*step)
                            assert enclosed(derivative,result[b][order]['gradient'][k]),(n,ustr,b,order,k,str(derivative),result[b][order]['gradient'][k])
                            gradients+=1
        if n%50==0:print('audited points',n+1,flush=True)
    receipt=dict(precision_digits=70,parameter_points=len(points),temperature_derivative_value_checks=checked,
                 gradient_checks=gradients,all_passed=True,seed=12102045,
                 method='Direct matrix exponential, temperature finite differences at step 1e-14, parameter/time finite differences at step 1e-12',
                 scope='Independent numerical cross-check of analytic first/second temperature derivatives and their interval gradients; not a box certificate')
    (ROOT/'temperature_reduction_independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
