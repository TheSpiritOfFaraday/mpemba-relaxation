"""Exact rational disk/half-plane inverse slow-mode feasibility.

No floating comparisons are used. A feasible return includes a strictly
positive reversible generator, an exact simple slow eigenvalue, and a positive
rational lower bound for the gap to both faster modes.
Completeness for the original closed triplet inequalities requires slack=0.
Positive slack strengthens those inequalities; rejection then excludes only
that strengthened problem, while a successful construction remains valid.
"""
import itertools,json
from pathlib import Path
import sympy as s

PAIRS=list(itertools.combinations(range(4),2))
TRIPLES=list(itertools.combinations(range(4),3))

def require(ok,msg):
    if not ok:raise ValueError(msg)

def exact_realize(q_input,w_input,slack=s.Rational(0)):
    slack=s.Rational(slack)
    require(slack>=0,'nonnegative rational slack')
    def result(**data):
        return dict(data,slack=str(slack),
            feasibility_domain='original triplet inequalities' if slack==0 else 'strengthened half-planes',
            infeasibility_excludes_original=(not data['feasible'] and slack==0),
            scope='Complete exact rational feasibility for the original premise' if slack==0 else
                'Exact strengthened feasibility; failure does not exclude the original premise')
    q=s.Matrix(list(map(s.Rational,q_input)));q/=sum(q)
    require(all(q[i]>q[i+1]>0 for i in range(3)),'strict ordered bath weights')
    w=s.Matrix(list(map(s.Rational,w_input)));w-=s.ones(4,1)*(q.T*w)[0]
    require(any(w),'nonconstant mode')
    D=s.diag(*q);null=(q.T.col_join((D*w).T)).nullspace()
    require(len(null)==2,'fast space dimension')
    U=s.Matrix.hstack(*null)
    bases=[s.eye(2)/2,s.diag(1,-1)/2,s.Matrix([[0,1],[1,0]])/2]
    down=s.Matrix([[-q[i]*(U[i,:]*B*U[j,:].T)[0] for B in bases] for i,j in PAIRS])
    rows=[]
    for i,j,k in TRIPLES:
        t=q[j]/q[i];u=q[k]/q[j]
        A=down[PAIRS.index((i,k)),:];B=down[PAIRS.index((j,k)),:];C=down[PAIRS.index((i,j)),:]
        rows.extend([A*(1+t*u/2)-B*(1+u/2)-C*(1-t),A*t-B])
    half=[]
    for row in rows:
        c,a,b=map(s.factor,list(row));c+=slack*max(abs(c),abs(a),abs(b))
        if a==b==0:
            if c>0:return result(feasible=False,reason='inconsistent constant half-plane')
        else:half.append((c,a,b))
    candidates={(s.Rational(0),s.Rational(0))}
    for c,a,b in half:candidates.add((-c*a/(a*a+b*b),-c*b/(a*a+b*b)))
    for (c,a,b),(f,d,e) in itertools.combinations(half,2):
        det=a*e-b*d
        if det:candidates.add(((b*f-c*e)/det,(c*d-a*f)/det))
    valid=[p for p in candidates if all(c+a*p[0]+b*p[1]<=0 for c,a,b in half)]
    if not valid:return result(feasible=False,reason='empty half-plane intersection',candidate_count=len(candidates))
    xy=min(valid,key=lambda p:p[0]**2+p[1]**2);radius=s.factor(xy[0]**2+xy[1]**2)
    if radius>=1:return result(feasible=False,reason='closest point outside open disk',min_radius_squared=str(radius),candidate_count=len(candidates))
    B=bases[0]+xy[0]*bases[1]+xy[1]*bases[2]
    require(B.det()>0 and s.trace(B)==1,'positive fast matrix')
    gram=U.T*D*U;gap=s.factor((B.det()*gram.det())/s.trace(B*gram))
    require(gap>0,'positive fast spectral lower bound')
    core=U*B*U.T;rho=max([s.Rational(0)]+[core[i,j] for i,j in PAIRS])+gap/100
    M=-D*core+rho*(q*s.ones(1,4)-s.eye(4))
    require(s.ones(1,4)*M==s.zeros(1,4),'column conservation')
    require(M*D==D*M.T,'detailed balance')
    require(M.T*w==-rho*w,'exact slow eigenpair')
    require(all(M[i,j]>0 for i in range(4) for j in range(4) if i!=j),'positive rates')
    evaluations=[]
    for i,j,k in TRIPLES:
        A=M[i,k];B0=M[j,k];C=M[i,j];t=q[j]/q[i];u=q[k]/q[j]
        F=A*(1+t*u/2)-B0*(1+u/2)-C*(1-t);G=A*t-B0;kappa=(A+B0)/(C*(1+t))
        P=(kappa>=1 and F>0) or (kappa<=1 and G>0)
        require(not P and F<=0 and G<=0,'original triplet predicate')
        evaluations.append(dict(triplet=[i+1,j+1,k+1],kappa=str(kappa),F=str(F),G=str(G),predicate=False))
    return result(feasible=True,verified=True,q=list(map(str,q)),w=list(map(str,w)),
        fast_basis=[[str(v) for v in U.row(i)] for i in range(4)],
        B=[[str(v) for v in B.row(i)] for i in range(2)],coordinates=list(map(str,xy)),
        min_radius_squared=str(radius),candidate_count=len(candidates),
        slow_eigenvalue=str(-rho),fast_gap_lower_bound=str(gap),
        generator=[[str(v) for v in M.row(i)] for i in range(4)],
        downward_rates=[str(M[i,j]) for i,j in PAIRS],triplets=evaluations)

if __name__=='__main__':
    q=['1','199/200','2/5','1/400'];w=['1','-759/1000','-48957/80000','-1/250']
    out=exact_realize(q,w,slack=s.Rational(1,10**8))
    require(out['feasible'],'known companion should be realizable')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['generator','triplets','downward_rates']},indent=2))
