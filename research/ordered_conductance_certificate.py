"""Exact weighted-contrast identities and primary-example compatibility.

The dimension-independent proof is in ordered_conductance_theorem.tex.
Symbolic checks through n=7 audit its formula; they do not replace that proof.
"""
import itertools,json,time
from pathlib import Path
import sympy as s
start=time.time()
def require(ok,msg):
    if not ok:raise ValueError(msg)
records=[]
for n in range(2,8):
    q=s.symbols('q0:'+str(n),positive=True);D=s.diag(*q)
    kk={ij:s.Symbol('k%d%d'%ij) for ij in itertools.combinations(range(n),2)}
    def k(i,j):return kk[tuple(sorted((i,j)))]
    L=s.zeros(n)
    for i,j in kk:L[i,j]=q[j]*k(i,j);L[j,i]=q[i]*k(i,j)
    for i in range(n):L[i,i]=-sum(L[i,j] for j in range(n) if i!=j)
    hs=[];etas=[]
    for m in range(n-1):
        Q=sum(q[:m+1]);h=s.Matrix([1 if i<=m else -Q/q[m+1] if i==m+1 else 0 for i in range(n)])
        hs.append(h);etas.append(Q*sum(q[:m+2])/q[m+1])
        require(s.cancel((h.T*D*s.ones(n,1))[0])==0,'zero mean')
        require(s.cancel((h.T*D*h)[0]-etas[-1])==0,'contrast norm')
    off=[]
    for m,l in itertools.combinations(range(n-1),2):
        require(s.cancel((hs[m].T*D*hs[l])[0])==0,'orthogonality')
        lhs=(hs[m].T*D*L*hs[l])[0]
        rhs=sum(q[:l+2])*sum(q[i]*(k(m+1,l+1)-k(i,l+1)) for i in range(m+1))
        rhs+=sum(q[j]*sum(q[i]*(k(m+1,j)-k(i,j)) for i in range(m+1)) for j in range(l+2,n))
        require(s.cancel(lhs-rhs)==0,'off diagonal identity')
        # Replace each column by its base and nonnegative successive increments.
        repl={k(i,j):s.Symbol('base%d'%j)+sum(s.Symbol('d%d_%d'%(r,j)) for r in range(1,i+1)) for i,j in kk}
        poly=s.Poly(s.expand(rhs.subs(repl)),*sorted(s.expand(rhs.subs(repl)).free_symbols,key=str))
        require(all(co>0 for co in poly.coeffs()),'positive polynomial coefficients')
        off.append(dict(indices=[m+1,l+1],positive_terms=len(poly.terms())))
    records.append(dict(n=n,off_diagonal_checks=off))

# Exact primary counterexample rates: all G<0, so this theorem permits its
# ordinary cooling reversal but excludes both strong and heating effects.
q=s.Matrix([1,s.Rational(199,200),s.Rational(2,5),s.Rational(1,400)]);q/=sum(q)
rates=list(map(s.Rational,['1','28','249/10','1799/8','4983/200','251/25']))
kk={(i,j):r/q[i] for (i,j),r in zip(itertools.combinations(range(4),2),rates)}
g=[]
for i,j,k0 in itertools.combinations(range(4),3):
    diff=kk[j,k0]-kk[i,k0];require(diff>0,'primary strict G ordering');g.append(str(diff))
out=dict(verified=True,symbolic_dimensions=list(range(2,8)),checks=records,
    primary_reduced_conductance_differences=g,elapsed_seconds=time.time()-start,
    scope='Exact identities and positive-coefficient decomposition; all-n proof is separate')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
