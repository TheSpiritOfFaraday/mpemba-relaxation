"""Exact rational enclosures of finite-time KL and TV crossing.
Uniformization with rigorously bounded Poisson tail; algebraic initial bounds
from the primary certificate; log series with rigorous remainder. All acceptance
arithmetic is fractions.Fraction, with directed rational rounding to 50 decimals.
"""
from fractions import Fraction as F
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1];r=json.loads((ROOT/'research/adversarial_certificate.json').read_text());SCALE=10**50

def rnd(iv):
 a,b=iv;return F(a.numerator*SCALE//a.denominator,SCALE),F(-((-b.numerator*SCALE)//b.denominator),SCALE)
def pt(x):return (F(x),F(x))
def add(a,b):return rnd((a[0]+b[0],a[1]+b[1]))
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[x*y for x in a for y in b];return rnd((min(v),max(v)))
def div(a,b):
 assert b[0]*b[1]>0;return mul(a,(1/b[1],1/b[0]))
def isum(xs):
 v=pt(0)
 for x in xs:v=add(v,x)
 return v
q=list(map(F,r['q']));Q=sum(q);pi=[z/Q for z in q];M=[[F(z) for z in row] for row in r['generator']];nu=F(253)
P=[[F(i==j)+M[i][j]/nu for j in range(4)] for i in range(4)]
assert all(z>=0 for row in P for z in row);assert all(sum(P[i][j] for i in range(4))==1 for j in range(4))
starts=[]
for key,n in [('hot_root_intervals',1000),('cold_root_intervals',100)]:
 z=[tuple(map(F,a)) for a in r[key]]
 for qi,(l,u) in zip(q,z):assert 0<=l and l**n<=qi<=u**n
 norm=isum(z);starts.append([div(a,norm) for a in z])

def evolve(t,start,N=250):
 x=nu*t;assert x<N+2
 # exp(x) enclosed via scalar Taylor series with positive remainder.
 term=F(1);exp_sum=F(1)
 for k in range(1,N+1):term*=x/k;exp_sum+=term
 nextterm=term*x/(N+1);tail=nextterm/(1-x/F(N+2))
 em=(1/(exp_sum+tail),1/exp_sum)
 state=start[:];tot=start[:];weight=pt(1)
 for k in range(1,N+1):
  state=[isum([mul(pt(P[i][j]),state[j]) for j in range(4)]) for i in range(4)]
  weight=mul(weight,pt(x/k));tot=[add(tot[i],mul(weight,state[i])) for i in range(4)]
 # P^k p has each entry in [0,1]; positive omitted sum bounded by scalar tail.
 return [mul(add(v,(F(0),tail)),em) for v in tot]

def logpoint(a,N=60):
 assert a>0;z=(a-1)/(a+1);az=abs(z);assert az<1
 power=pt(z);z2=pt(z*z);total=pt(0)
 for k in range(N):
  total=add(total,mul(power,pt(F(2,2*k+1))));power=mul(power,z2)
 rem=2*az**(2*N+1)/((2*N+1)*(1-az*az))
 return add(total,(-rem,rem))
def logiv(a):return(logpoint(a[0])[0],logpoint(a[1])[1])
def absiv(a):return(max(F(0),a[0],-a[1]),max(abs(a[0]),abs(a[1])))
def distances(p):
 tv=mul(pt(F(1,2)),isum([absiv(sub(z,pt(pii))) for z,pii in zip(p,pi)]))
 kl=isum([mul(z,logiv(div(z,pt(pii)))) for z,pii in zip(p,pi)])
 return {'TV':tv,'KL':kl}

out={'method':'exact rational uniformization and log Taylor bounds','decimal_rounding':50,'uniformization_rate':'253','poisson_terms':251,'log_terms':60,'times':{}}
for time in ['1/10','11/100','14/100','16/100']:
 t=F(time);ds=[distances(evolve(t,p)) for p in starts];delta={k:sub(ds[0][k],ds[1][k]) for k in ['TV','KL']}
 out['times'][time]={'hot':{k:list(map(str,v)) for k,v in ds[0].items()},'cold':{k:list(map(str,v)) for k,v in ds[1].items()},'hot_minus_cold':{k:list(map(str,v)) for k,v in delta.items()}}
 print(time,{k:[float(x) for x in v] for k,v in delta.items()},flush=True)
assert F(out['times']['1/10']['hot_minus_cold']['TV'][0])>0
assert F(out['times']['11/100']['hot_minus_cold']['TV'][1])<0
assert F(out['times']['14/100']['hot_minus_cold']['KL'][0])>0
assert F(out['times']['16/100']['hot_minus_cold']['KL'][1])<0
out['verified']='TV crossing in (0.10,0.11), KL crossing in (0.14,0.16) by continuity; no uniqueness claim'
(ROOT/'research/finite_time_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('ALL EXACT FINITE-TIME CROSSING CHECKS PASSED')
