"""Original-field soft branches, hard type sums and collapsed spin contractions.
Independent soft extraction: Bessel split of Hankel modes then helicity Wick
contractions. Compared with the 0F1 field endpoints, without any basis rotation.
"""
from itertools import product
from pathlib import Path
import mpmath as mp
import numpy as np
from aoki_majorana import I, EPSILON_LOWER as EPS, J, slash, contraction
from aoki_weyl import endpoints, soft_branch, hard_block, scalar_coefficients
mp.mp.dps = 40
lines = []


def check(name, residual, tolerance=5e-12):
    assert residual < tolerance, (name, residual)
    line=f'PASS {name}: residual {float(residual):.3e}'
    lines.append(line);print(line,flush=True)


def bessel_modes(mu, x):
    # Split each Hankel function into its two J branches before any contraction.
    z=1j*mu;nf=mp.sqrt(mp.pi*x)/2*mp.exp(mp.pi*mu/2)
    ng=mp.sqrt(mp.pi*x)/2*mp.exp(-mp.pi*mu/2)
    nu=.5-z
    fp=nf*mp.besselj(-nu,x)/(1j*mp.sin(mp.pi*nu))
    fm=-nf*mp.exp(-1j*mp.pi*nu)*mp.besselj(nu,x)/(1j*mp.sin(mp.pi*nu))
    nu=.5+z
    gm=ng*mp.besselj(-nu,x)/(1j*mp.sin(mp.pi*nu))
    gp=-ng*mp.exp(-1j*mp.pi*nu)*mp.besselj(nu,x)/(1j*mp.sin(mp.pi*nu))
    return {1:(fp,gp),-1:(fm,gm)}


def bessel_soft(mu,k,t1,t2,branch):
    p=np.linalg.norm(k)
    f1,g1=bessel_modes(mu,-p*t1)[branch]
    f2,g2=bessel_modes(mu,-p*t2)[-branch]
    out=np.zeros((4,4),complex)
    phase=mp.exp(3j*mp.pi/4)/mp.sqrt(2)
    for h in (1,-1):
        u1,v1=map(complex,(phase*(f1-h*g1),phase*(f1+h*g1)))
        u2,v2=map(complex,(phase*(f2-h*g2),phase*(f2+h*g2)))
        proj=(I+h*slash(k/p))/2
        out[:2,:2]-=u1*v2.conjugate()*proj@EPS
        out[:2,2:]+=u1*u2.conjugate()*proj
        out[2:,:2]-=v1*v2.conjugate()*EPS@proj@EPS
        out[2:,2:]+=v1*u2.conjugate()*EPS@proj
    return out


def scalar_hankel(mu,p,t1,t2,a,b):
    def modes(t):
        x=-p*t
        return [complex(mp.sqrt(mp.pi*x)/2*mp.exp(s*mp.pi*mu/2)*
                        mp.hankel1(.5-s*1j*mu,x)) for s in (1,-1)]
    f1,g1=modes(t1);f2,g2=modes(t2)
    gt=np.array([g1*f2.conjugate(),f1*f2.conjugate(),g1*g2.conjugate()])
    lt=np.array([f1.conjugate()*g2,-g1.conjugate()*g2,-f1.conjugate()*f2])
    if (a,b)==(-1,1):return gt
    if (a,b)==(1,-1):return lt
    th=float(t1>t2)+.5*float(t1==t2)
    return th*gt+(1-th)*lt if a==1 else (1-th)*gt+th*lt

rng=np.random.default_rng(20260930)
for mu in (.2,.73,2.,5.):
    errors=[]
    for scale in (.04,.4,1.3):
        k=scale*np.array([.3,-.4,.7])
        for branch in (1,-1):
            d=soft_branch(mu,k,-.37,-1.41,branch)
            expected=bessel_soft(mu,k,-.37,-1.41,branch)
            errors.append(np.linalg.norm(d-expected)/np.linalg.norm(expected))
    check(f'both soft branches: 0F1 vs Bessel/helicity contractions, mu={mu}',max(errors))

for mu in (.2,.73,2.):
    radial_errors=[];cut_errors=[];zero_errors=[];mixed_errors=[]
    for trial in range(2):
        p=rng.normal(size=3);n=p/np.linalg.norm(p)
        q=.1*rng.normal(size=3);K=.1*rng.normal(size=3)
        ts=-rng.uniform(.15,2.,4)
        for a,b in product((1,-1),repeat=2):
            d=contraction(p,ts[0],ts[1],a,b,mu)
            wgf,wff,wgg=scalar_hankel(mu,np.linalg.norm(p),ts[0],ts[1],a,b)
            radial_errors.append(np.linalg.norm(np.array(scalar_coefficients(d,n))-[wgf,wff,wgg]))
            left0=np.vstack((I,EPS));right0=np.hstack((EPS,I))
            leftm=np.vstack((I,-EPS));rightm=np.hstack((EPS,-I))
            v0=-.5*right0@J@d@J@left0
            dneg=contraction(-p,ts[0],ts[1],a,b,mu)
            vneg=-.5*right0@J@dneg@J@left0
            zero_errors.extend([np.linalg.norm(v0+wgf*slash(n)),np.linalg.norm(v0+vneg)])
            mixed_errors.extend([np.linalg.norm(-.5*rightm@J@d@J@left0-wff*I),
                                 np.linalg.norm(-.5*right0@J@d@J@leftm+wgg*I)])
        p2=rng.normal(size=3)
        for aa in product((1,-1),repeat=4):
            d1=contraction(p,ts[0],ts[1],aa[0],aa[1],mu)
            d2=contraction(p2,ts[2],ts[3],aa[2],aa[3],mu)
            for sq,sK in product((1,-1),repeat=2):
                direct=np.trace(J@d1@J@soft_branch(mu,K,ts[1],ts[2],sK)@
                                J@d2@J@soft_branch(mu,q,ts[3],ts[0],sq))
                left=hard_block(d1,mu,q,K,ts[0],ts[1],(sq,sK))
                right=hard_block(d2,mu,K,q,ts[2],ts[3],(sK,sq))
                c=complex(2**(-1-2j*mu)*mp.gamma(.5-1j*mu)**2/mp.pi)
                cp={1:c,-1:-c.conjugate()}
                pref=cp[sq]*cp[sK]*(np.dot(q,q)*ts[3]*ts[0])**(sq*1j*mu)
                pref*= (np.dot(K,K)*ts[1]*ts[2])**(sK*1j*mu)
                fact=pref*np.trace(left@slash(K/np.linalg.norm(K))@right@slash(q/np.linalg.norm(q)))
                cut_errors.append(abs(direct-fact)/max(abs(direct),1e-15))
    check(f'four-type scalar coefficients vs Hankel, all SK, mu={mu}',max(radial_errors))
    check(f'full type-cycle vs factorized cut, all SK and soft assignments, mu={mu}',max(cut_errors),2e-10)
    check(f'zeroth hard block and Bose oddness, mu={mu}',max(zero_errors))
    check(f'mixed scalar hard blocks from original field types, mu={mu}',max(mixed_errors))
    # Local first descendants at generic hard direction, all SK orderings.
    errors=[]
    q0=np.array([.31,-.22,.14]);K0=q0+np.array([0.,0.,1.])
    for scale in (.004,.002,.001):
        err=0.
        for a,b in product((1,-1),repeat=2):
            d=contraction(p,-.37,-1.41,a,b,mu)
            wgf,wff,wgg=scalar_hankel(mu,np.linalg.norm(p),-.37,-1.41,a,b)
            q=scale*q0;K=scale*K0
            exact=hard_block(d,mu,q,K,-.37,-1.41)
            approx=-wgf*slash(n)+1j*(-.37)/(1+2j*mu)*wff*slash(q)
            approx-=1j*(-1.41)/(1+2j*mu)*wgg*slash(K)
            err=max(err,np.linalg.norm(exact-approx))
        errors.append(err)
    assert errors[1]<.26*errors[0] and errors[2]<.26*errors[1],errors
    lines.append(f'PASS first descendants converge quadratically, mu={mu}: '+', '.join(f'{e:.3e}' for e in errors))
    print(lines[-1],flush=True)

lines.append('Scope: finite-time spin/field-type identities and local soft expansion; not full BD box quadrature.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
