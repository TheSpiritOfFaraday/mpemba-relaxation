"""Adversarial checks of the exact inverse constructor and its output.

The output generator is checked by standard-library rational LDL inertia,
separately from the constructor's 2x2 fast-matrix argument.
"""
from fractions import Fraction as F
import json,time
from pathlib import Path
import sympy as s
from exact_mode_realizability import exact_realize
from standalone_primary_check import row_generator,positive_inertia,require
start=time.time();checks=[]
q=['1','199/200','2/5','1/400'];w=['1','-759/1000','-48957/80000','-1/250']
base=exact_realize(q,w)
require(base['feasible'],'known feasible mode')
for label,ww in [('original',w),('scaled',[str(-7*s.Rational(v)) for v in w]),
                 ('translated',[str(s.Rational(v)+31) for v in w])]:
 out=exact_realize(q,ww)
 require(out['feasible'],label+' affine mode invariance')
 qq=list(map(F,out['q']));rates=list(map(F,out['downward_rates']))
 L=row_generator(qq,rates);lam=F(out['slow_eigenvalue']);gap=F(out['fast_gap_lower_bound'])
 upper=lam+min(-lam/2,gap/4);lower=lam-gap/2
 require(positive_inertia(L,qq,upper)[0]==1,label+' independently counted slow upper bracket')
 require(positive_inertia(L,qq,lower)[0]==2,label+' independently counted slow lower bracket')
 ww=list(map(F,out['w']))
 require(all(sum(L[i][j]*ww[j] for j in range(4))==lam*ww[i] for i in range(4)),label+' exact row eigenpair')
 checks.append(dict(test=label,passed=True))
# This Gibbs overlap has three prescribed positive roots, including the bath;
# the all-G theorem excludes its two nonbath strong roots under triplet failure.
strong=exact_realize(['1','1/2','1/4','1/8'],['-3/32','11/16','-3/2','1'])
require(not strong['feasible'],'strong-effect direction excluded')
checks.append(dict(test='strong direction excluded',passed=True,receipt=strong))
# First contrast violates F whenever it is the simple slow direction here.
first=exact_realize(['1','1/4','1/16','1/64'],['1','-4','0','0'])
require(not first['feasible'],'first contrast excluded by full triplet premise')
checks.append(dict(test='first contrast excluded',passed=True,receipt=first))
strengthened=exact_realize(q,w,slack=2)
require(not strengthened['feasible'],'large slack excludes this strengthened problem')
require(not strengthened['infeasibility_excludes_original'],'positive slack rejection is not original infeasibility')
require(first['infeasibility_excludes_original'],'zero slack rejection addresses original premise')
checks.append(dict(test='strengthened rejection scope',passed=True,receipt=strengthened))
try:exact_realize(q,w,slack=-1)
except ValueError:checks.append(dict(test='negative slack rejected',passed=True))
else:raise ValueError('negative slack should be rejected')
for label,qq,ww in [('equal energies',['1','1','1/2','1/3'],w),
                    ('zero bath weight',['1','1/2','1/3','0'],w),
                    ('constant mode',q,['2']*4)]:
 try:exact_realize(qq,ww)
 except ValueError:checks.append(dict(test=label,passed=True))
 else:raise ValueError(label+' should be rejected')
out=dict(verified=True,checks=checks,elapsed_seconds=time.time()-start,
 scope='Adversarial inputs, affine invariance, and independent rational inertia of accepted constructions')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
