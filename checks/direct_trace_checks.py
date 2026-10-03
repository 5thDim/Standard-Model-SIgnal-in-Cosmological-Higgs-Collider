"""Direct traces with exact external factors; no dummy-pair average.

Run: PYTHONDONTWRITEBYTECODE=1 python3 -u checks/direct_trace_checks.py
Each time and SK label stays attached to its physical external vertex.
"""
from itertools import product
import mpmath as mp
import numpy as np
from aoki_majorana import contraction, slash
from aoki_weyl import hard_block, scalar_coefficients
mp.mp.dps=25

def check(name,error,tol):
    assert error<tol,(name,error,tol)
    print(f"PASS {name}: residual {float(error):.3e}",flush=True)

def coeff(mu,k,t1,t2,a,b):
    n=np.array([0.,0.,1.])
    return scalar_coefficients(contraction(k*n,t1,t2,a,b,mu),n)[:2]

def raw_vectors(mu,k,n,q,K,t1,t2,a,b,right=False):
    k0,k1=coeff(mu,k,t1,t2,a,b)
    reverse=coeff(mu,k,t2,t1,b,a)[1]
    h=1e-4*k
    v=[coeff(mu,k+j*h,t1,t2,a,b)[0] for j in (-2,-1,1,2)]
    radial=k*(v[0]-8*v[1]+8*v[2]-v[3])/(12*h)-k0
    ell=(q+K)/2
    hard=k0*ell+radial*n*(n@ell)
    A=t1*k1;B=t2*reverse
    qin,qout=(K,q) if right else (q,K)
    w12=hard-1j*k/(1+2j*mu)*(A*qin+B*qout)
    w21=hard-1j*k/(1+2j*mu)*(B*qin+A*qout)
    p1=k0-1j*k/(1+2j*mu)*(A+B)
    compact=p1*ell+radial*n*(n@ell)
    return w12,w21,compact

def physical_block(mu,k,n,q,K,t1,t2,a,b,reverse=False,right=False):
    ell=(q+K)/2
    qin,qout=(K,q) if right else (q,K)
    if reverse:
        # Vertex 2 is now adjacent to qin; labels remain physical labels.
        t1,t2,a,b=t2,t1,b,a
        n=-n
    d=contraction(k*n+ell,t1,t2,a,b,mu)
    return hard_block(d,mu,qin,qout,t1,t2)

def E(k,t,a):
    return (1-1j*a*k*t)*np.exp(1j*a*k*t)

n=np.array([np.sin(.7),0.,np.cos(.7)])
m=np.array([np.sin(1.1)*np.cos(.8),np.sin(1.1)*np.sin(.8),np.cos(1.1)])
q0=np.array([.31,-.22,.14]);K0=q0+np.array([0.,0.,1.])
for mu in (.2,.73,2.):
    firsterrors=[];reverseerrors=[]
    for a,b in product((-1,1),repeat=2):
        k=1.3;t1,t2=-.43,-1.27;h=1e-4
        w12,w21,compact=raw_vectors(mu,k,n,q0,K0,t1,t2,a,b)
        for reverse,w in ((False,w12),(True,w21)):
            derivative=(physical_block(mu,k,n,h*q0,h*K0,t1,t2,a,b,reverse)
                        -physical_block(mu,k,n,-h*q0,-h*K0,t1,t2,a,b,reverse))/(2*h)
            firsterrors.append(np.linalg.norm(derivative+slash(w)/k)/np.linalg.norm(derivative))
        reverseerrors.append(abs(coeff(mu,k,t1,t2,a,b)[0]-coeff(mu,k,t2,t1,b,a)[0]))
        assert np.linalg.norm(w12+w21-2*compact)<1e-12
    check(f"separate 12 and 21 first-order blocks, all SK, mu={mu}",max(firsterrors),2e-7)
    check(f"zeroth-order reversal identity, all SK, mu={mu}",max(reverseerrors),1e-13)

errors=[];weightederrors=[]
mu=.73;q=.017*q0;K=.017*K0;ell=(q+K)/2;ks=K-q
kL,kR=1.3,.9
times=(-.43,-1.27,-.61,-1.49)
energies=(np.linalg.norm(kL*n+ks/2),np.linalg.norm(-kL*n+ks/2),
          np.linalg.norm(kR*m-ks/2),np.linalg.norm(-kR*m-ks/2))
assert energies[0]!=energies[1] and energies[2]!=energies[3]
for a,b,c,d in product((-1,1),repeat=4):
    left=raw_vectors(mu,kL,n,q,K,times[0],times[1],a,b)
    right=raw_vectors(mu,kR,m,q,K,times[2],times[3],c,d,True)
    trace=0j
    scalar=0j
    for wl,wr in product(left[:2],right[:2]):
        trace+=np.trace((-slash(wl)/kL)@slash(K/np.linalg.norm(K))
                        @(-slash(wr)/kR)@slash(q/np.linalg.norm(q)))
        scalar+=2/(kL*kR*np.linalg.norm(q)*np.linalg.norm(K))*(
            2*(wl@ell)*(wr@ell)-.5*(wl@ks)*(wr@ks)
            -(wl@wr)*(ell@ell-ks@ks/4))
    vl,vr=left[2],right[2]
    compact=8/(kL*kR*np.linalg.norm(q)*np.linalg.norm(K))*(
        2*(vl@ell)*(vr@ell)-.5*(vl@ks)*(vr@ks)
        -(vl@vr)*(ell@ell-ks@ks/4))
    errors.extend([abs(trace-scalar),abs(trace-compact)])
    ext=np.prod([E(k,t,s) for k,t,s in zip(energies,times,(a,b,c,d))])
    weightederrors.append(abs(ext*trace-ext*compact))
check("four individual traces vs scalar expression, all 16 SK",max(errors),1e-12)
check("same equality with common exact unequal Higgs factors",max(weightederrors),1e-12)

# Check first surviving full soft-cut trace without averaging any time integrand.
for mu in (.2,.73,2.):
    errors=[]
    for scale in (.008,.004,.002):
        q=scale*q0;K=scale*K0;ell=(q+K)/2;ks=K-q
        a,b,c,d=1,-1,-1,1
        vl=sum(physical_block(mu,kL,n,q,K,times[0],times[1],a,b,rev) for rev in (False,True))
        vr=sum(physical_block(mu,kR,m,q,K,times[2],times[3],c,d,rev,True) for rev in (False,True))
        raw=np.trace(vl@slash(K/np.linalg.norm(K))@vr@slash(q/np.linalg.norm(q)))
        wl=raw_vectors(mu,kL,n,q,K,times[0],times[1],a,b)[2]
        wr=raw_vectors(mu,kR,m,q,K,times[2],times[3],c,d,True)[2]
        leading=np.trace((-2*slash(wl)/kL)@slash(K/np.linalg.norm(K))
                         @(-2*slash(wr)/kR)@slash(q/np.linalg.norm(q)))
        errors.append(abs(raw-leading)/abs(leading))
    # At fixed, unsymmetrized times the omitted second-order hard-end
    # corrections give an O(soft) relative remainder. No time average is used.
    assert errors[1]<.6*errors[0] and errors[2]<.6*errors[1],errors
    print(f"PASS direct four-diagram soft-cut expansion mu={mu}: "+
          ", ".join(f"{e:.3e}" for e in errors),flush=True)
print("Scope: fixed-time original-field contractions and four actual diagrams, with exact external factors; no dummy-pair averaging, no full BD-box integration.",flush=True)
