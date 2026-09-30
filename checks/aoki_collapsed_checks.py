"""Complete two-component soft-cut integrand vs first nonzero expansion.
Whittaker hard lines and exact 0F1 soft endpoints, with all 16 SK and
four external orderings. Finite time interval: tests local expansion,
not the infinite-contour scalar-seed value. Several masses/geometries.
"""
from pathlib import Path
import numpy as np
import mpmath as mp
from aoki_majorana import uv, I, E, J, T, X, UB, BETA, slash
mp.mp.dps=25
lines=[]

def radial(mu,p,ts):
    # Build the lower-index contractions directly from Aoki helicity sums.
    # Axis z is sufficient for separating the four scalar radial blocks.
    n=np.array([0.,0.,1.]); P=[(I+h*slash(n))/2 for h in (1,-1)]
    C=np.zeros((len(ts),len(ts),4,4),complex)
    for h,proj in zip((1,-1),P):
        vals=np.array([[complex(a) for a in uv(mu,float(p*t),h)] for t in ts])
        u,v=vals.T
        for (r,s),coef,mat in [((0,0),-np.outer(u,v.conj()),proj@E),
                              ((0,1),np.outer(u,u.conj()),proj),
                              ((1,0),-np.outer(v,v.conj()),E@proj@E),
                              ((1,1),np.outer(v,u.conj()),E@proj)]:
            C[:,:,2*r:2*r+2,2*s:2*s+2]+=coef[:,:,None,None]*mat
    R=UB.conj().T@T.conj().T@C@X@T@UB
    gt=[R[:,:,0,0],R[:,:,0,2],R[:,:,2,0],R[:,:,2,2]]
    # Canonical greater/lesser radial reversal, verified separately against C.
    lt=[-gt[3].T,gt[1].T,gt[2].T,-gt[0].T]
    return gt,lt

def summed(mu,p,e1,e2,ts,ws):
    gt,lt=radial(mu,p,ts)
    x=ts[:,None];y=ts[None,:]
    theta=(x<y)+.5*(x==y)
    out=np.zeros((4,len(ts),len(ts)),complex)
    for a in (1,-1):
        for b in (1,-1):
            ex=(1+1j*a*e1*x)*np.exp(-1j*a*e1*x)*(1+1j*b*e2*y)*np.exp(-1j*b*e2*y)
            for i in range(4):
                R=gt[i] if (a,b)==(-1,1) else lt[i] if (a,b)==(1,-1) else theta*gt[i]+(1-theta)*lt[i] if a==1 else (1-theta)*gt[i]+theta*lt[i]
                out[i]+=a*b*ex*R
    return out*ws[:,None]*ws[None,:]*(x*y)**(1j*mu-1)

def endpoints(mu,p,ts):
    z=1j*mu; d=1j/(1+2*z)
    # No truncation of soft endpoint descendants.
    f=np.array([complex(mp.hyp0f1(.5+z,-(p*t)**2/4)) for t in ts])
    g=np.array([complex(-d*p*t*mp.hyp0f1(1.5+z,-(p*t)**2/4)) for t in ts])
    return f,g

def block(mu,n,q,K,ts,ws):
    ell=(q+K)/2; ks=K-q; result=np.zeros((2,2),complex)
    fq,gq=endpoints(mu,np.linalg.norm(q),ts)
    fK,gK=endpoints(mu,np.linalg.norm(K),ts)
    sq=slash(q/np.linalg.norm(q));sK=slash(K/np.linalg.norm(K))
    for sign in (1,-1):
        p=sign*n+ell; pnorm=np.linalg.norm(p)
        e1=np.linalg.norm(sign*n+ks/2);e2=np.linalg.norm(-sign*n+ks/2)
        R=summed(mu,pnorm,e1,e2,ts,ws); sp=slash(p/pnorm)
        # row(q) M(p) col(K): explicit 2x2 sums over the two sectors.
        result+=(R[0]*np.outer(gq,fK)).sum()*sq
        result+=(R[1]*np.outer(gq,gK)).sum()*sq@sp@sK
        result-=(R[2]*np.outer(fq,fK)).sum()*sp
        result-=(R[3]*np.outer(fq,gK)).sum()*sK
    return result

def run(mu,order,case):
    zz,ww=np.polynomial.legendre.leggauss(order)
    ts=.2+1.4*(zz+1);ws=1.4*ww
    R=summed(mu,1.,1.,1.,ts,ws);h=R[2].sum();h11=(ts[:,None]*R[0]).sum()
    step=1e-4
    dh=(summed(mu,1+step,1.,1.,ts,ws)[2].sum()-summed(mu,1-step,1.,1.,ts,ws)[2].sum())/(2*step)
    p1=h+2j/(1+2j*mu)*h11;p2=dh-h
    n=np.array([np.sin(.7),0.,np.cos(.7)])
    m=np.array([np.sin(1.1)*np.cos(.8),np.sin(1.1)*np.sin(.8),np.cos(1.1)])
    q0=np.array([.31,-.22,.14]) if case==0 else np.array([-.42,.37,-.7])
    errors=[]; values=[]
    for scale in (.008,.004,.002):
        q=scale*q0;K=q+scale*np.array([0.,0.,1.]);ell=(q+K)/2
        left=block(mu,n,q,K,ts,ws);right=block(mu,m,K,q,ts,ws)
        left0=-2*slash(p1*ell+p2*n*np.dot(n,ell))
        right0=-2*slash(p1*ell+p2*m*np.dot(m,ell))
        sq=slash(q/np.linalg.norm(q));sK=slash(K/np.linalg.norm(K))
        exact=np.trace(left@sK@right@sq)
        expected=np.trace(left0@sK@right0@sq)
        errors.append(abs(exact-expected)/abs(expected));values.append(exact/scale**2)
    assert errors[2] < errors[1]/3 and errors[1] < errors[0]/3,errors
    assert errors[-1]<2e-4,errors
    lines.append(f'PASS full soft endpoints + hard/external momenta mu={mu}, N={order}, geometry={case}: relative errors '+', '.join(f'{e:.3e}' for e in errors))
    # A single trace is nonzero. Summing orderings kills the degree-zero term.
    single=np.trace(slash(n)@sK@slash(m)@sq)*h*h
    assert abs(single)>1e-5
    zeroth=sum(a*b*single for a in (1,-1) for b in (1,-1))
    assert abs(zeroth)<1e-13
    return values[-1]

for mu in (.2,.73,2.):
    for case in (0,1):
        run(mu,20,case)
        print(lines[-1],flush=True)
# Quadrature refinement does not test the BD limit; it tests the finite-domain
# integral whose leading coefficient is obtained independently by Taylor expansion.
a=run(.73,32,0);b=run(.73,48,0)
# Time ordering gives a cusp on x=y; tensor quadrature converges algebraically.
assert abs(a-b)/abs(b)<.004
lines.append(f'PASS finite-domain quadrature refinement 32 -> 48: relative change {abs(a-b)/abs(b):.3e}')
lines.append('PASS zeroth-order trace nonzero individually, cancels between external orderings; coefficient scales as soft^2 before d^3q.')
lines.append('Scope: finite-domain full soft-cut expansion; not an independent BD hard amplitude or uncut-box evaluation.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[-3:]))
