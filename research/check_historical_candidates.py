"""Check retained discovery claims without treating finite searches as proofs.

Exact arithmetic checks the auxiliary nonmonotone example and two backward
designs. Direct 100-digit evolution rechecks the broad heating search's best
point. Terminal logs authenticate recorded counters only, not excluded points.
"""
from pathlib import Path
import hashlib
import itertools
import json
import sympy as s
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
checks = []
inputs = set()


def require(ok, label):
    if not ok:
        raise ValueError(label)
    checks.append(label)


def read(name):
    inputs.add('research/' + name)
    return json.loads((ROOT / 'research' / name).read_text())


def generator(q, rates):
    n = len(q)
    M = s.zeros(n)
    for (i, j), a in zip(itertools.combinations(range(n), 2), rates):
        M[i, j] = a
        M[j, i] = a * q[j] / q[i]
    for j in range(n):
        M[j, j] = -sum(M[i, j] for i in range(n) if i != j)
    return M


def triplets(M, q):
    rows = []
    for i, j, k in itertools.combinations(range(len(q)), 3):
        A, B, C = M[i, k], M[j, k], M[i, j]
        t, u = q[j]/q[i], q[k]/q[j]
        H = A+B-C*(1+t)
        G = A*t-B
        F = A*(1+t*u/2)-B*(1+u/2)-C*(1-t)
        rows.append((H, G, F))
    return rows


def main():
    R = s.Rational
    q = s.Matrix([1, R(3, 4), R(1, 2), R(1, 4)])
    M = generator(q, [R(1, 4), R(1, 16), R(3, 4), R(5, 4), 6, 9])
    expected = s.Matrix([[-R(13,32),R(1,4),R(1,16),R(3,4)],
                         [R(3,16),-R(37,12),R(5,4),6],
                         [R(1,32),R(5,6),-R(93,16),9],
                         [R(3,16),2,R(9,2),-R(63,4)]])
    require(M == expected, 'Auxiliary printed generator exactly reconstructed')
    require(s.ones(1,4)*M == s.zeros(1,4), 'Auxiliary column conservation')
    require(M*s.diag(*q) == s.diag(*q)*M.T, 'Auxiliary detailed balance')
    rows = triplets(M, q)
    require(rows == [(R(7,8),-R(77,64),-R(317,192)),
                     (R(101,16),-R(87,16),-R(199,32)),
                     (R(309,32),-R(69,8),-R(167,16)),
                     (R(155,12),-5,-R(14,3))], 'All twelve printed H/G/F values')
    require(all(h>0 and g<0 and f<0 for h,g,f in rows), 'Strict upper-branch failures')
    x = s.Symbol('x')
    P = 1536*x**3+38480*x**2+175584*x+100755
    require(s.expand(M.charpoly(x).as_expr()-x*P/1536) == 0, 'Printed characteristic polynomial')
    for a,b in [(-20,-19),(-6,-5),(-R(67,100),-R(66,100))]:
        require(P.subs(x,a)*P.subs(x,b)<0, 'Sign-changing cubic bracket '+str((a,b)))
    # Three disjoint brackets exhaust a cubic's roots and imply simplicity.
    w = (x*s.eye(4)-M.T).adjugate()[:,0]
    expected_w = s.Matrix([x**3+R(1183,48)*x**2+R(10031,96)*x+R(6717,256),
                           x*x/4+R(1333,192)*x+R(6717,256),
                           x*x/16+R(467,96)*x+R(6717,256),
                           3*x*x/4+R(559,64)*x+R(6717,256)])
    require(s.simplify(w-expected_w)==s.zeros(4,1), 'Printed cofactor eigenfunction')
    diffs = [-x*(192*x*x+4684*x+18729)/192,
             -x*(12*x+133)/64, x*(132*x+743)/192]
    require(all(s.expand(w[i+1]-w[i]-diffs[i])==0 for i in range(3)), 'Adjacent difference identities')
    lo = -R(67,100)
    require(4684*lo+18729>0 and 12*lo+133>0 and 132*lo+743>0,
            'Exact interval signs are positive, positive, negative; no monotone orientation')

    q = s.Matrix([1,R(199,200),R(2,5),R(1,400)])
    w = s.Matrix([1,-R(759,1000),-R(48957,80000),-R(1,250)])
    lam = -R(4966671,250000)
    a,b = s.symbols('a b')
    designs = []
    for variant, rates in [('rejected',[1,28,R(249,10),R(1799,8),a,b]),
                            ('accepted',[1,28,R(249,10),a,R(4983,200),b])]:
        G = generator(q,rates)
        sol = s.solve(list(G.T*w-lam*w),(a,b))
        G = G.subs(sol)
        require(G.T*w == lam*w, variant+' exact prescribed eigenpair')
        fg = triplets(G,q)
        positive = [(i+1,str(f),str(g)) for i,(_,g,f) in enumerate(fg) if f>0 or g>0]
        require(bool(positive)==(variant=='rejected'), variant+' original triplet verdict')
        if variant=='accepted':
            certificate=read('adversarial_rational_mode_certificate.json')
            require(G==s.Matrix([[R(v) for v in row] for row in certificate['generator']]),
                    'Accepted design matches the certified rational companion')
        designs.append({'design':variant,'positive_triplet_tests':positive})

    # The apparent negative arithmetic margin has a stronger exact rejection:
    # subtracting the minimum reset leaves exactly two connected components.
    saved=read('general_n_arithmetic_search_result.json')['best']
    qn=s.Matrix([R(str(v)) for v in saved['q']]); qn/=sum(qn)
    require(all(qn[i]>qn[i+1]>0 for i in range(len(qn)-1)),
            'Borderline candidate has positive strictly ordered exact bath weights')
    kn=list(map(lambda v:R(str(v)),saved['reduced_rates']))
    n=len(qn); pairn=list(itertools.combinations(range(n),2))
    require(n==6 and min(kn)==-1, 'Borderline arithmetic candidate has exact reset level')
    core=s.zeros(n); adjacency=[set() for _ in range(n)]
    for (i,j),k in zip(pairn,kn):
        c=k+1
        require(c>=0, 'Nonnegative reset-subtracted edge '+str((i+1,j+1)))
        core[i,j]=qn[i]*c; core[j,i]=qn[j]*c
        if c>0: adjacency[i].add(j); adjacency[j].add(i)
    for j in range(n):core[j,j]=-sum(core[i,j] for i in range(n) if i!=j)
    remaining=set(range(n));components=[]
    while remaining:
        todo=[min(remaining)];component=set()
        while todo:
            i=todo.pop()
            if i in component:continue
            component.add(i);todo.extend(adjacency[i]-component)
        remaining-=component;components.append(sorted(component))
    require(components==[[0,1,2,3],[4,5]], 'Exactly two connected reset-subtracted components')
    reset=qn*s.ones(1,n)-s.eye(n)
    require(all(g<=0 and f<=0 for _,g,f in triplets(core+reset,qn)),
            'Borderline candidate satisfies every original F/G inequality exactly')
    mode=s.Matrix([sum(qn[4:])]*4+[-sum(qn[:4])]*2)
    require(core.T*mode==s.zeros(n,1) and (qn.T*mode)[0]==0,
            'Exact component-constant mode is bath-orthogonal and core-stationary')
    require((core+reset).T*mode==-mode and mode[0]+mode[1]-2*mode[2]==0,
            'Exact slow eigenpair has zero arithmetic margin, not a negative violation')
    # Reversible Dirichlet positivity makes the connected-component core's
    # other eigenvalues strictly negative. Two core zero modes leave exactly
    # one bath-orthogonal zero mode; adding reset proves it is simple and slow.

    status=read('transient_heating_status.json')
    mp.mp.dps=100
    p=status['best'];E=list(map(lambda z:mp.mpf(str(z)),p['E']))
    qm=mp.matrix([mp.exp(-e) for e in E]);qm/=sum(qm)
    K=mp.matrix([[mp.mpf(str(v)) for v in row] for row in p['K']])
    G=mp.matrix(len(E))
    for j in range(len(E)):
        for i in range(len(E)):
            if i!=j:G[i,j]=qm[i]*K[i,j]
        G[j,j]=-sum(G[i,j] for i in range(len(E)) if i!=j)
    beta=mp.mpf(str(p['beta']));t=mp.mpf(str(p['time']))
    initial=mp.matrix([mp.exp(-beta*e) for e in E]);initial/=sum(initial)
    mean=sum(initial[i]*E[i] for i in range(len(E)))
    derivative_initial=mp.matrix([-initial[i]*(E[i]-mean) for i in range(len(E))])
    U=mp.expm(t*G);pt=U*initial;dp=U*derivative_initial
    terms=[dp[i]*mp.log(pt[i]/qm[i]) for i in range(len(E))]
    derivative=sum(terms);score=derivative/sum(abs(v) for v in terms)
    require(derivative>0 and score>0, 'Broad heating best point remains positive at100 digits')
    require(abs(score-mp.mpf(str(p['score'])))<mp.mpf('1e-8'), 'Broad heating saved score numerical correspondence')
    require(status['evaluations']%(200*16*24)==0, 'Broad search counter counts full attempted batches')
    terminal={}
    for stem in ['transient_heating','transient_heating_optimized']:
        name=stem+'_search.log' if stem=='transient_heating' else stem+'.log'
        inputs.add('research/'+name)
        line=(ROOT/'research'/name).read_text().splitlines()[-1]
        require(line.startswith('FINISHED '), stem+' terminal log exists')
        saved=read(stem+'_status.json')
        require(json.loads(line[len('FINISHED '):])==saved, stem+' terminal log matches status')
        require(saved['status']=='finished' and saved['hits']==0 and saved['best_score']>0,
                stem+' retained positive best score and zero accepted hits')
        terminal[stem]={'evaluations':saved['evaluations'],'elapsed_seconds':saved['elapsed_seconds']}
    optimized=read('transient_heating_optimized_status.json')
    require(optimized['evaluations']==(optimized['run']+1)*501*12*11,
            'Optimizer call count matches1004 completed fixed-size runs, including penalties')
    inputs.add('research/theorem_notes.md')
    result={'passed':True,'count':len(checks),'checks':checks,'backward_designs':designs,
            'broad_best_point':{'precision':100,'derivative':str(derivative),'normalized_score':str(score)},
            'exact_roundoff_rejection':{'components':components,'slow_eigenvalue':'-1','arithmetic_margin':'0','analytic_dependency':'Reversible connected-component Dirichlet form plus reset shift proves simple slow ordering.'},
            'terminal_records':terminal,'scope':'Exact auxiliary/design and two-component roundoff-rejection algebra; numerical saved-point check; historical terminal-counter correspondence. No all-domain heating conclusion and no reauthentication of excluded samples.',
            'input_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(inputs)},
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'research/detailed_audit/historical_candidate_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
