"""Independent exact arithmetic for the focused paper's worked example.

Standard library only. No imported project code, saved certificate, numerical
eigenvector, floating-point acceptance, or assert statements are used.
The 4x4 characteristic polynomial is reconstructed from principal minors.
Positive rational nth-root intervals use integer bisection, not float guesses.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
checks = []

def require(ok, label):
    if not ok:
        raise ValueError(label)
    checks.append(label)

def determinant(a):
    n = len(a)
    answer = F(0)
    for p in permutations(range(n)):
        term = F((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n)))
        for i in range(n):
            term *= a[i][p[i]]
        answer += term
    return answer

def add(a,b): return a[0]+b[0], a[1]+b[1]
def mul(a,b):
    terms = [x*y for x in a for y in b]
    return min(terms), max(terms)
def divide(a,b):
    require(b[0] > 0, 'interval denominator positive')
    return mul(a, (1/b[1],1/b[0]))
def squared(a):
    return (F(0) if a[0]<=0<=a[1] else min(x*x for x in a), max(x*x for x in a))

def nth_root(q,n):
    scale = 10**15
    low, high = 0, scale
    target = q.numerator * scale**n
    while high-low>1:
        mid = (low+high)//2
        if mid**n*q.denominator <= target:
            low = mid
        else:
            high = mid
    result = F(low,scale), F(high,scale)
    require(result[0]**n<=q<=result[1]**n, f'integer root enclosure n={n}, q={q}')
    return result

def show(interval):
    return {'exact':[str(x) for x in interval], 'decimal_display_only':[float(x) for x in interval]}

def main():
    q = list(map(F,['1','199/200','2/5','1/400']))
    rates = list(map(F,['1','28','249/10','330154675693/1470375000','4983/200','1527327066/151990625']))
    w = list(map(F,['1','-759/1000','-48957/80000','-1/250']))
    lam = -F(4966671,250000)
    M = [[F(0) for j in range(4)] for i in range(4)]
    for (i,j),rate in zip(combinations(range(4),2),rates):
        M[i][j] = rate
        M[j][i] = rate*q[j]/q[i]
    for j in range(4):
        M[j][j] = -sum(M[i][j] for i in range(4) if i!=j)
        require(sum(M[i][j] for i in range(4))==0, f'column {j+1} conserves mass')
    for i,j in combinations(range(4),2):
        require(M[i][j]>0 and M[j][i]>0, f'edge {i+1}{j+1} positive in both directions')
        require(M[i][j]*q[j]==M[j][i]*q[i], f'edge {i+1}{j+1} balances')
    Q=sum(q); pi=[x/Q for x in q]
    require(pi==list(map(F,['400/959','398/959','160/959','1/959'])), 'printed equilibrium')
    require(all(sum(M[i][j]*pi[j] for j in range(4))==0 for i in range(4)), 'stationarity')
    triplets=[]
    bounds=[('-230','-196'),('-1/50','-13/100'),('-19/10','-2/25'),('-118','-3/100')]
    for (i,j,k),(fb,gb) in zip(combinations(range(4),3),bounds):
        A,B,C=M[i][k],M[j][k],M[i][j]
        t,u=q[j]/q[i],q[k]/q[j]
        kap=(A+B)/(C*(1+t)); H=A+B-C*(1+t)
        ff=A*(1+t*u/2)-B*(1+u/2)-C*(1-t); gg=A*t-B
        require(ff==(2/(1+t)+u/2)*gg+(1-t)/(1+t)*H, 'branch-reduction identity')
        require(not ((kap>=1 and ff>0) or (kap<=1 and gg>0)), 'original inclusive-branch predicate fails')
        require(ff<F(fb) and gg<F(gb), 'printed strict triplet bounds')
        triplets.append({'states':f'{i+1}{j+1}{k+1}', 'kappa':str(kap),'F':str(ff),'G':str(gg)})
    require(sum(qi*wi for qi,wi in zip(q,w))==0, 'bath-orthogonal mode')
    require(all(sum(w[i]*M[i][j] for i in range(4))==lam*w[j] for j in range(4)), 'left eigenpair')
    N=sum(q[i]*w[i]**2 for i in range(4)); v=[q[i]*w[i]/N for i in range(4)]
    require(N==F(27567998009,16000000000),'printed mode normalization')
    require(sum(w[i]*v[i] for i in range(4))==1,'paired normalization')
    require(all(sum(M[i][j]*v[j] for j in range(4))==lam*v[i] for i in range(4)), 'right eigenpair')
    # det(xI-M): coefficient of x^(4-k) is (-1)^k times sum of principal k-minors.
    coeff=[F(1)]
    for k in range(1,5):
        coeff.append((-1)**k*sum(determinant([[M[i][j] for j in subset] for i in subset])
                                      for subset in combinations(range(4),k)))
    trace=sum(M[i][i] for i in range(4)); s=trace-lam
    trace2=sum(M[i][j]*M[j][i] for i in range(4) for j in range(4))
    p=(s*s-(trace2-lam*lam))/2
    require(coeff==[1,-s-lam,p+lam*s,-lam*p,0], 'full characteristic factor x(x-lambda)(x^2-s*x+p)')
    poly=lambda x: x*x-s*x+p
    for lo,hi in [(-337,-336),(-60,-59)]:
        require(poly(lo)*poly(hi)<0, f'one fast root in ({lo},{hi})')
    require(-20<lam<-19, 'simple slow eigenvalue in (-20,-19); gap greater than 39')
    def preparation(n):
        z=[nth_root(qi,n) for qi in q]
        Z=(sum(x[0] for x in z),sum(x[1] for x in z))
        probs=[divide(x,Z) for x in z]
        numerator=(F(0),F(0))
        for wi,zi in zip(w,z): numerator=add(numerator,mul((wi,wi),zi))
        overlap=divide(numerator,Z)
        distance=(F(0),F(0))
        for i in range(4):
            d=(probs[i][0]-pi[i],probs[i][1]-pi[i])
            distance=add(distance,divide(squared(d),(q[i],q[i])))
        return overlap,distance,probs
    hot,Dh,ph=preparation(1000); cold,Dc,pc=preparation(100)
    require(F('-.093755177882')<hot[0]<=hot[1]<F('-.093755177881'),'printed hot bracket')
    require(F('-.093857287853')<cold[0]<=cold[1]<F('-.093857287852'),'printed cold bracket')
    require(cold[1]<hot[0] and hot[1]<0,'nonzero reversed overlap ordering')
    require(Dc[1]<24<Dh[0], 'initial hotter-farther ordering')
    ratio=divide((-hot[1],-hot[0]),(-cold[1],-cold[0]))
    advantage=(1-ratio[1],1-ratio[0])
    result={'passed':True,'scope':'Exact focused-witness arithmetic; eventual dynamics use the separately checked spectral proof.',
            'checks':checks,'q':list(map(str,q)),'rates':list(map(str,rates)),
            'generator':[[str(x) for x in row] for row in M], 'triplets':triplets,
            'lambda':str(lam),'fast_quadratic_coefficients':list(map(str,[1,-s,p])),
            'hot_overlap':show(hot),'cold_overlap':show(cold),
            'hot_initial_Dq':show(Dh),'cold_initial_Dq':show(Dc),'amplitude_advantage':show(advantage),
            'hot_probabilities':[show(x) for x in ph],'cold_probabilities':[show(x) for x in pc],
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'paper_sha256':hashlib.sha256((ROOT/'focused_paper.tex').read_bytes()).hexdigest()}
    out=ROOT/'research/detailed_audit/freshman_arithmetic.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS: {len(checks)} exact checks; {out}')
    print('Amplitude advantage (display only):',show(advantage)['decimal_display_only'])

if __name__=='__main__': main()
