import numpy as np
from itertools import combinations
q=np.array([1,199/200,2/5,1/400]); a=np.array([-24,3,-.1,200,.04,.04]);pairs=list(combinations(range(4),2));E=-np.log(q);p=q/q.sum()
a+=25*q[[i for i,j in pairs]]
M=np.zeros((4,4))
for z,(i,j) in enumerate(pairs):M[i,j]=a[z];M[j,i]=a[z]*q[j]/q[i]
M-=np.diag(M.sum(axis=0)); S=M*np.sqrt(p)[None,:]/np.sqrt(p)[:,None];l,V=np.linalg.eigh(S);w=V[:,-2]/np.sqrt(p)
print('q',q,'E',E,'a',a,'M',M,'l',l,'w',w)
for i,j,k in combinations(range(4),3):
 A=M[i,k];B=M[j,k];C=M[i,j];t=q[j]/q[i];u=q[k]/q[j]
 print((i,j,k),(A+B)/(C*(1+t)),A*(1+t*u/2)-B*(1+u/2)-C*(1-t),A*t-B)
for b in [0,.0001,.001,.002,.005,.01,.02,.05,.1]:
 pp=q**b;pp/=pp.sum();ov=w@(pp-p);der=-w@(pp*(E-pp@E));print(b,ov,der,ov*der)
