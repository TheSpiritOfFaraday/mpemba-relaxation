"""Two-dimensional convex feasibility for a prescribed four-state slow mode.
Floating discovery implementation. Exact rational implementation/proof separate.
"""
import itertools
import numpy as np
from cone_search import constraints,pairs

def realize(q,w,tol=1e-12,slack=0.):
    q=np.asarray(q,float);q=q/sum(q);w=np.asarray(w,float);w=w-w@q
    if max(abs(w))==0:raise ValueError('Nonconstant mode required')
    w=w/max(abs(w))
    _,_,vh=np.linalg.svd(np.array([q,q*w]),full_matrices=True);F=vh[2:].T
    BB=np.array([np.eye(2)/2,np.diag([1.,-1.])/2,np.array([[0.,1.],[1.,0.]])/2])
    aa=np.array([[-q[i]*(F[i]@B@F[j]) for i,j in pairs] for B in BB])
    rows=constraints(q)@aa.T
    sizes=np.max(abs(rows),axis=1)
    rows=rows[sizes>0]/sizes[sizes>0,None]
    c=rows[:,0]+slack;A=rows[:,1:]
    norm=np.linalg.norm(A,axis=1);keep=norm>1e-15
    if np.any(c[~keep]>tol):return dict(feasible=False,min_radius_squared=float('inf'))
    A=A[keep]/norm[keep,None];c=c[keep]/norm[keep]
    candidates=[np.zeros(2)]+[-ci*ai for ai,ci in zip(A,c)]
    for i,j in itertools.combinations(range(len(c)),2):
        mm=A[[i,j]]
        if abs(np.linalg.det(mm))>1e-12:candidates.append(np.linalg.solve(mm,-c[[i,j]]))
    valid=[p for p in candidates if np.max(A@p+c,initial=-np.inf)<=tol]
    if not valid:return dict(feasible=False,min_radius_squared=float('inf'))
    z=min(valid,key=lambda p:p@p);radius=float(z@z)
    if radius>=1-1e-8:return dict(feasible=False,min_radius_squared=radius,coordinates=z.tolist())
    B=BB[0]+z[0]*BB[1]+z[1]*BB[2]
    # Reset adds no F/G residual and restores all rates to strict positivity.
    rr=aa[0]+z[0]*aa[1]+z[1]*aa[2]
    reset=q[[i for i,j in pairs]];rho=max(0,float(np.max(-rr/reset)))+1
    rates=rr+rho*reset
    return dict(feasible=True,min_radius_squared=radius,coordinates=z.tolist(),rates=rates.tolist(),reset=rho,F=F.tolist(),B=B.tolist(),w=w.tolist())

def contrast_basis(q):
    q=np.asarray(q,float);q=q/sum(q);h=np.zeros((4,3))
    for k in range(3):
        Q=sum(q[:k+1]);eta=Q*sum(q[:k+2])/q[k+1]
        h[:k+1,k]=1/np.sqrt(eta);h[k+1,k]=-Q/q[k+1]/np.sqrt(eta)
    return h

if __name__=='__main__':
    import json
    q=[1,199/200,2/5,1/400];w=[1,-.75893729814493,-.61211930644096,-.00386630776426]
    print(json.dumps(realize(q,w),indent=2))
