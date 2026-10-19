"""Discovery-only multistart search of rank <=3 self-Chollet ratios."""
import itertools, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from numba import njit

PERMS = np.array(list(itertools.permutations(range(5))), dtype=np.int64)

@njit(cache=True)
def per_and_grad(C):
    p=0.0
    g=np.zeros((5,5))
    for perm in PERMS:
        a=np.empty(5)
        for i in range(5): a[i]=C[i,perm[i]]
        t=1.0
        for i in range(5): t*=a[i]
        p+=t
        for i in range(5):
            t=1.0
            for j in range(5):
                if j != i: t*=a[j]
            g[i,perm[i]]+=t
    return p,g

@njit(cache=True)
def objective(x):
    V=x.reshape((5,3))
    C=V@V.T
    p,dp=per_and_grad(C)
    q,dq=per_and_grad(C*C)
    r=q/(p*p)
    dC=2*C*dq/(p*p)-2*r*dp/p
    dV=(dC+dC.T)@V
    return -r,-dV.ravel()

def norm(V):
    return V/np.linalg.norm(V,axis=1)[:,None]

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--starts',type=int,default=500)
    ap.add_argument('--seed',type=int,default=20261014)
    args=ap.parse_args()
    rng=np.random.default_rng(args.seed)
    out=Path(__file__).parent
    best=0.; records=[]; t0=time.time()
    for k in range(args.starts):
        V=norm(rng.normal(size=(5,3)))
        opt=minimize(objective,V.ravel(),method='BFGS',jac=True,
                     options={'gtol':1e-10,'maxiter':1500})
        W=norm(opt.x.reshape(5,3)); value=-objective(W.ravel())[0]
        records.append({'start':k,'ratio':float(value),'success':bool(opt.success),'nit':int(opt.nit)})
        if value>best:
            best=value
            np.savez(out/'best.npz',V=W,C=W@W.T,ratio=value)
            print(k,'best',best,'elapsed',time.time()-t0,flush=True)
        if value>13/48+1e-8:
            print('NUMERICAL COUNTEREXAMPLE',W,W@W.T,flush=True)
            break
    (out/'search_log.json').write_text(json.dumps({'seed':args.seed,'elapsed':time.time()-t0,'records':records},indent=2))
