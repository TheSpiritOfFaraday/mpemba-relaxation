"""Exact finite-temperature witness accompanying the continuous family proof.

Only integer arithmetic and fractions are used for thermal weights.  The
generator, gap, and triplet signs are established by quarter_limit_certificate.
No previously generated certificate is read as input here.
"""
from fractions import Fraction as F
from pathlib import Path
import json,time
from standalone_primary_check import approximate_overlap,require
start=time.time()
e=F(1,10**50)
q=(F(1),1-e/3,F(3,5),e*e)
w=(F(1),F(-4,5),-F(1,3)-4*e/9+5*e**3/3,-e)
require(sum(qi*wi for qi,wi in zip(q,w))==0,'exact bath orthogonality')
hot=approximate_overlap(q,w,10000)
cold=approximate_overlap(q,w,50)
require(hot[1]<0 and cold[1]<0,'negative overlaps')
gain=(1-(-hot[0])/(-cold[1]),1-(-hot[1])/(-cold[0]))
require(gain[0]>F(22,100),'amplitude advantage exceeds 22 percent')
# The analytic perturbation lemma proves a uniform overlap change <=1000 e^3
# after adding e^4 H, where all reduced conductances of H are <=100.
err=1000*e**3
strict_hot=(hot[0]-err,hot[1]+err)
strict_cold=(cold[0]-err,cold[1]+err)
strict_gain=1-(-strict_hot[0])/(-strict_cold[1])
require(strict_cold[1]<0 and strict_hot[1]<0,'strict-family overlap signs')
require(strict_gain>F(22,100),'strict-triplet advantage exceeds 22 percent')
out=dict(verified=True,epsilon=str(e),q=list(map(str,q)),w=list(map(str,w)),
 beta_hot='1/10000',beta_cold='1/50',beta_bath='1',
 overlap_hot=list(map(str,hot)),overlap_cold=list(map(str,cold)),
 advantage_enclosure=list(map(str,gain)),
 displayed_advantage_percent=[100*float(v) for v in gain],
 perturbation_overlap_error_bound=str(err),strict_advantage_lower_bound=str(strict_gain),
 elapsed_seconds=time.time()-start,
 scope='Thermal interval certificate; generator proof and perturbation bound are separate analytic dependencies')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['q','w','overlap_hot','overlap_cold','advantage_enclosure','strict_advantage_lower_bound']},indent=2))
