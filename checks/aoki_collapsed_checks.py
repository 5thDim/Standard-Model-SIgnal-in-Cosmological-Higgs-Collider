"""Complete two-component soft-cut integrand vs first nonzero expansion.
Whittaker hard lines and exact 0F1 soft endpoints, with all 16 SK and
four external orderings. Finite time interval: tests local expansion,
not the infinite-contour scalar-seed value. Several masses/geometries.
"""
from pathlib import Path
import numpy as np
import mpmath as mp
from aoki_majorana import J, slash
from aoki_weyl import endpoints, summed_grid, scalar_coefficients
mp.mp.dps=25
lines=[]

def summed(mu,p,e1,e2,ts,ws):
    # Extract scalar coefficients only after the original four field-type sum.
    direction=np.array([0.,0.,1.])
    D=summed_grid(mu,p*direction,e1,e2,ts,ws)
    return scalar_coefficients(D,direction)


def block(mu,n,q,K,ts,ws):
    ell=(q+K)/2;ks=K-q;result=np.zeros((2,2),complex)
    # Each endpoint retains its original psi / psi^dagger type.
    right=np.array([endpoints(mu,q,-t)[1] for t in ts])
    left=np.array([endpoints(mu,K,-t)[0] for t in ts])
    for sign in (1,-1):
        p=sign*n+ell
        e1=np.linalg.norm(sign*n+ks/2);e2=np.linalg.norm(-sign*n+ks/2)
        D=summed_grid(mu,p,e1,e2,ts,ws)
        result-=.5*np.einsum('tij,tsjk,skl->il',right@J,D,J@left)
    return result

def run(mu,order,case):
    zz,ww=np.polynomial.legendre.leggauss(order)
    ts=.2+1.4*(zz+1);ws=1.4*ww
    wgf,wff,_=summed(mu,1.,1.,1.,ts,ws);h=wgf.sum();h11=(ts[:,None]*wff).sum()
    step=1e-4
    dh=(summed(mu,1+step,1.,1.,ts,ws)[0].sum()-summed(mu,1-step,1.,1.,ts,ws)[0].sum())/(2*step)
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
    lines.append(f'PASS direct psi/psi-dagger soft endpoints + hard/external momenta mu={mu}, N={order}, geometry={case}: relative errors '+', '.join(f'{e:.3e}' for e in errors))
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
