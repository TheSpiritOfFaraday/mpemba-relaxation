"""Exact rational structural checks and high-precision thermal screening.
Uses the existing certified quarter-limit family; no historical files modified.
The numerical conclusions are not interval certificates or hardware validation.
"""
import ast,itertools,json
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
mp.mp.dps=160
ROOT=Path(__file__).resolve().parents[1]
source=json.loads((ROOT/'quarter_limit_certificate.json').read_text())
pairs=list(itertools.combinations(range(4),2))
def parse(s,e):
 def go(n):
  if isinstance(n,ast.Constant):return F(n.value)
  if isinstance(n,ast.Name) and n.id=='e':return e
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -go(n.operand)
  if isinstance(n,ast.BinOp):
   a,b=go(n.left),go(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
   if isinstance(n.op,ast.Pow):
    assert b.denominator==1
    return a**int(b)
  raise ValueError(ast.dump(n))
 return go(ast.parse(s,mode='eval').body)
def m(x):return mp.mpf(x.numerator)/x.denominator
def calc(e,delta):
 q=[F(1),1-e/3,F(3,5),e*e];pi=[v/sum(q) for v in q]
 M=[[parse(s,e) for s in row] for row in source['generator']]
 H=[[F(0) for _ in range(4)] for _ in range(4)]
 for (i,j),k in zip(pairs,[100,10,1,11,2,3]):H[i][j]=pi[i]*k;H[j][i]=pi[j]*k
 for j in range(4):H[j][j]=-sum(H[i][j] for i in range(4) if i!=j)
 M=[[M[i][j]+delta*H[i][j] for j in range(4)] for i in range(4)]
 assert all(sum(M[i][j] for i in range(4))==0 for j in range(4))
 assert all(M[i][j]*pi[j]==M[j][i]*pi[i] for i,j in pairs)
 assert all(M[i][j]>0 for i in range(4) for j in range(4) if i!=j)
 residuals=[]
 for i,j,k in itertools.combinations(range(4),3):
  A,B,C=M[i][k],M[j][k],M[i][j];t=q[j]/q[i];u=q[k]/q[j]
  ff=A*(1+t*u/2)-B*(1+u/2)-C*(1-t);gg=A*t-B
  assert ff<=0 and gg<=0
  if delta:assert ff<0 and gg<0
  residuals.append(dict(triplet=[i+1,j+1,k+1],F=str(ff),G=str(gg)))
 p=mp.matrix(list(map(m,pi)));sqrtp=[mp.sqrt(x) for x in p]
 S=mp.matrix([[m(M[i][j])*sqrtp[j]/sqrtp[i] for j in range(4)] for i in range(4)])
 lam,V=mp.eigsy(S);w=mp.matrix([V[i,2]/sqrtp[i] for i in range(4)])
 def prep(beta):
  v=[mp.power(m(x),beta) for x in q];return mp.matrix([x/sum(v) for x in v])
 hot,cold=prep(mp.mpf('0.0001')),prep(mp.mpf('0.09'))
 bh,bc=mp.fdot(w,hot),mp.fdot(w,cold)
 gain=1-abs(bh/bc)
 def deviation(v,t):
  coeff=V.T*mp.matrix([(v[i]-p[i])/sqrtp[i] for i in range(4)])
  z=V*mp.matrix([coeff[i]*mp.exp(lam[i]*t) if i<3 else 0 for i in range(4)])
  return mp.matrix([z[i]*sqrtp[i] for i in range(4)])
 def kl(d):
  # Subtract the exactly zero total linear term to avoid cancellation.
  return sum((p[i]+d[i])*mp.log1p(d[i]/p[i])-d[i] for i in range(4))
 def distances(tau):
  t=mp.mpf(tau)/(-lam[2]);return kl(deviation(hot,t)),kl(deviation(cold,t))
 times=[]
 for tau in [0,.5,1,2,3,5,10,20,30,40,50,60,70,80,100]:
  dh,dc=distances(tau)
  times.append(dict(slow_time_units=tau,hot_KL=str(dh),cold_KL=str(dc),KL_advantage=str(1-dh/dc)))
 # A bracketed crossing is descriptive numerical evidence, not a certificate
 # of uniqueness or of the earliest crossing.
 low,high=mp.mpf(0),mp.mpf(100)
 if distances(high)[0]<distances(high)[1]:
  for _ in range(150):
   mid=(low+high)/2;dh,dc=distances(mid)
   if dh>dc:low=mid
   else:high=mid
  crossing=dict(slow_time_units=str((low+high)/2),KL=str(distances((low+high)/2)[1]))
 else:crossing=None
 return dict(epsilon=str(e),perturbation=str(delta),q=list(map(str,q)),generator=[[str(x) for x in r] for r in M],triplet_residuals=residuals,eigenvalues=[str(x) for x in lam],hot_beta='1/10000',cold_beta='9/100',bath_beta='1',amplitude_advantage=str(gain),capacitance_ratio=str(max(q)/min(q)),finite_time=times,numerical_crossing=crossing,structural_check='EXACT rational positivity, conservation, detailed balance and F/G signs',thermal_status='160-digit numerical evaluation, not interval-certified')
if __name__=='__main__':
 rows=[calc(F(1,1000),d) for d in [F(0),F(1,10**12),F(1,1000000),F(1,10000),F(1,1000)]]
 (Path(__file__).parent/'intermediate_results.json').write_text(json.dumps(rows,indent=2)+'\n')
 for r in rows:
  print('delta',r['perturbation'],'amplitude gain',float(r['amplitude_advantage']),'cap ratio',r['capacitance_ratio'])
  print('crossing',r['numerical_crossing'])
