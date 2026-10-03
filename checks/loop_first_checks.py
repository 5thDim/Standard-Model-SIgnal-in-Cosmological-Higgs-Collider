"""Loop-first Chapter 3 checks: no full BD-box quadrature is claimed.

Run: PYTHONDONTWRITEBYTECODE=1 python3 -u checks/loop_first_checks.py
Tests dummy-time symmetrization before integration, soft power counting,
the completed spin trace, and the loop coefficients with independent
left/right time kernels.
"""
from itertools import product
import mpmath as mp
import numpy as np
import sympy as sp
from aoki_majorana import contraction, slash, J
from aoki_weyl import hard_block, scalar_coefficients, soft_branch
import whittaker_scalar_seed as seeds

mp.mp.dps=25


def check(name, error, tol):
    assert error<tol, (name,error,tol)
    print(f"PASS {name}: residual {float(error):.3e}",flush=True)


def coefficients(mu,k,t1,t2,a,b):
    n=np.array([0.,0.,1.])
    return scalar_coefficients(contraction(k*n,t1,t2,a,b,mu),n)[:2]


def E(k,t,a):
    return (1-1j*a*k*t)*np.exp(1j*a*k*t)


def pi(mu,k,t1,t2,a,b):
    k0,k1=coefficients(mu,k,t1,t2,a,b)
    reverse=coefficients(mu,k,t2,t1,b,a)[1]
    h=1e-4*k
    # Five-point derivative at fixed times, independent of the seed bootstrap.
    vals=[coefficients(mu,k+j*h,t1,t2,a,b)[0] for j in (-2,-1,1,2)]
    dk=(vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*h)
    return k0-1j*k/(1+2j*mu)*(t1*k1+t2*reverse), k*dk-k0


def weighted_block(mu,k,n,q,K,t1,t2,a,b):
    """Unexpanded soft endpoints and hard line, two external orderings."""
    ell=(q+K)/2; ks=K-q
    result=np.zeros((2,2),complex)
    for eta in (-1,1):
        p=eta*k*n+ell
        e1=np.linalg.norm(eta*k*n+ks/2)
        e2=np.linalg.norm(-eta*k*n+ks/2)
        d=contraction(p,t1,t2,a,b,mu)
        result+=E(e1,t1,a)*E(e2,t2,b)*hard_block(d,mu,q,K,t1,t2)
    return result


n=np.array([np.sin(.7),0.,np.cos(.7)])
q0=np.array([.31,-.22,.14]); K0=q0+np.array([0.,0.,1.])
for mu in (.2,.73,2.):
    allerrors=[[],[],[]]
    for a,b in product((-1,1),repeat=2):
        k=1.3;t1,t2=-.43,-1.27
        p1,p2=pi(mu,k,t1,t2,a,b)
        ell0=(q0+K0)/2
        target=-2/k*slash(p1*ell0+p2*n*(n@ell0))*E(k,t1,a)*E(k,t2,b)
        assert np.linalg.norm(target)>1e-5
        for j,scale in enumerate((.008,.004,.002)):
            q,K=scale*q0,scale*K0
            # Symmetrize dummy time/SK pairs, not the value after integration.
            actual=(weighted_block(mu,k,n,q,K,t1,t2,a,b)
                    +weighted_block(mu,k,n,q,K,t2,t1,b,a))/2
            allerrors[j].append(np.linalg.norm(actual/scale-target)/np.linalg.norm(target))
    errors=list(map(max,allerrors))
    assert errors[1]<errors[0]/3 and errors[2]<errors[1]/3,errors
    check(f"pointwise dummy-symmetrized first-order block, all SK, mu={mu}",
          errors[-1],3e-5)
    print("  soft-halving relative errors: "+", ".join(f"{x:.3e}" for x in errors),flush=True)

# The three displayed terms immediately after Eq. (47), BEFORE dummy averaging.
def first_order_terms(mu,k,n,q,K,t1,t2,a,b):
    ell=(q+K)/2
    k0,k1=coefficients(mu,k,t1,t2,a,b)
    reverse=coefficients(mu,k,t2,t1,b,a)[1]
    _,p2=pi(mu,k,t1,t2,a,b)
    e1,e2=E(k,t1,a),E(k,t2,b)
    de1=k*t1*t1*np.exp(1j*a*k*t1)
    de2=k*t2*t2*np.exp(1j*b*k*t2)
    hard=-e1*e2/k*slash(k0*ell+p2*n*(n@ell))
    soft=1j*e1*e2/(1+2j*mu)*(t1*k1*slash(q)+t2*reverse*slash(K))
    external=-(n@(K-q))/2*(de1*e2-e1*de2)*k0*slash(n)
    return hard,soft,external


def single_ordering(mu,k,n,q,K,t1,t2,a,b):
    ell=(q+K)/2;ks=K-q
    p=k*n+ell
    e1=np.linalg.norm(k*n+ks/2)
    e2=np.linalg.norm(-k*n+ks/2)
    d=contraction(p,t1,t2,a,b,mu)
    return E(e1,t1,a)*E(e2,t2,b)*hard_block(d,mu,q,K,t1,t2)


raw_errors=[];average_errors=[];external_errors=[]
for mu in (.2,.73,2.):
    for a,b in product((-1,1),repeat=2):
        k=1.3;t1,t2=-.43,-1.27;step=1e-4
        terms=first_order_terms(mu,k,n,q0,K0,t1,t2,a,b)
        reverse_terms=first_order_terms(mu,k,n,q0,K0,t2,t1,b,a)
        derivative=(single_ordering(mu,k,n,step*q0,step*K0,t1,t2,a,b)
                    -single_ordering(mu,k,n,-step*q0,-step*K0,t1,t2,a,b))/(2*step)
        raw_errors.append(np.linalg.norm(sum(terms)-derivative)/np.linalg.norm(derivative))
        p1,p2=pi(mu,k,t1,t2,a,b);ell=(q0+K0)/2
        compact=-E(k,t1,a)*E(k,t2,b)/k*slash(p1*ell+p2*n*(n@ell))
        averaged=(sum(terms)+sum(reverse_terms))/2
        average_errors.append(np.linalg.norm(averaged-compact)/np.linalg.norm(compact))
        external_errors.append(np.linalg.norm(terms[2]+reverse_terms[2]))
check("three raw first-order sources vs unexpanded single-ordering derivative",
      max(raw_errors),2e-7)
check("dummy average of all three sources vs Pi1,Pi2",max(average_errors),2e-11)
check("external-mode first-order term cancels under dummy average",
      max(external_errors),1e-13)

# Cancellation requires external orderings, not the type sum of one diagram.
m=np.array([np.sin(1.1)*np.cos(.8),np.sin(1.1)*np.sin(.8),np.cos(1.1)])
single=np.trace(slash(n)@slash(K0/np.linalg.norm(K0))@slash(m)@slash(q0/np.linalg.norm(q0)))
assert abs(single)>1e-3
check("zeroth-order integrand: sum of four external orderings",
      abs(sum(a*b*single for a,b in product((-1,1),repeat=2))),1e-14)

# Full field-type sum of tr(DJDJDJDJ), without time or momentum integration.
mu=.73
q=.017*q0;K=.017*K0;ell=(q+K)/2
times=(-.43,-1.27,-.61,-1.49)
d23=soft_branch(mu,K,times[1],times[2])
d41=soft_branch(mu,q,times[3],times[0])
c=complex(mp.power(2,-1-2j*mu)*mp.gamma(.5-1j*mu)**2/mp.pi)
prefactor=c*c*(np.linalg.norm(q)*np.linalg.norm(K))**(2j*mu)*np.prod([(-t)**(1j*mu) for t in times])
errors=[]
for a,b,csk,d in product((-1,1),repeat=4):
    dl=contraction(1.3*n+ell,times[0],times[1],a,b,mu)
    dr=contraction(.9*m+ell,times[2],times[3],csk,d,mu)
    full=np.trace(dl@J@d23@J@dr@J@d41@J)
    vl=hard_block(dl,mu,q,K,times[0],times[1])
    vr=hard_block(dr,mu,K,q,times[2],times[3])
    factored=prefactor*np.trace(vl@slash(K/np.linalg.norm(K))@vr@slash(q/np.linalg.norm(q)))
    errors.append(abs(full-factored)/max(abs(full),1e-12))
check("full tr(DJDJDJDJ) vs soft factorization, all 16 SK",max(errors),1e-12)

# The scalar trace identity must hold for arbitrary unequal complex kernels.
rng=np.random.default_rng(1234)
err=[]
for _ in range(20):
    pl=rng.normal(size=2)+1j*rng.normal(size=2)
    pr=rng.normal(size=2)+1j*rng.normal(size=2)
    q=rng.normal(size=3);ks=rng.normal(size=3);K=q+ks;ell=(q+K)/2
    vl=pl[0]*ell+pl[1]*n*(n@ell)
    vr=pr[0]*ell+pr[1]*m*(m@ell)
    raw=np.trace((-2*slash(vl))@slash(K/np.linalg.norm(K))
                 @(-2*slash(vr))@slash(q/np.linalg.norm(q)))
    scalar=8/(np.linalg.norm(q)*np.linalg.norm(K))*(
        2*(vl@ell)*(vr@ell)-.5*(vl@ks)*(vr@ks)-(vl@vr)*(ell@ell-ks@ks/4))
    err.append(abs(raw-scalar)/max(1,abs(raw)))
check("complete trace with unequal complex left/right kernels",max(err),1e-12)

# Polynomial azimuthal moments, before any time integration.
u,v=sp.symbols("u v")
w=(v-u)/2
ll=(u+v)/2-sp.Rational(1,4)
rr=ll-w*w
polys={
    "M":ll*ll+ll/4-w*w/2,
    "P":(ll+sp.Rational(1,4))*rr/2,
    "Q":(ll+sp.Rational(1,4))*(w*w-rr/2)-w*w/2,
    "E":rr/2,
    "F":w*w-rr/2,
    "T":rr*rr/8,
}
polys["U"]=rr*w*w/2-polys["T"]
polys["V"]=w**4-6*polys["U"]-3*polys["T"]
polys["X"]=(ll-sp.Rational(1,4))*rr/2
polys["Y"]=(ll-sp.Rational(1,4))*(w*w-rr/2)
c1,c3=map(mp.mpf,map(str,(n[2],m[2])))
e=mp.mpf(str(n@m))
f0=lambda a,b:mp.gamma(a+b-mp.mpf("1.5"))*mp.gamma(mp.mpf("1.5")-a)*mp.gamma(mp.mpf("1.5")-b)/(mp.gamma(a)*mp.gamma(b)*mp.gamma(3-a-b))
maxerr=mp.mpf(0)
loop_values={}
for mu in map(mp.mpf,(".2",".73","2")):
    nu=mp.mpf(".5")-1j*mu
    base=f0(nu,nu)/(4*mp.pi)**mp.mpf("1.5")
    moments={}
    for name,poly in polys.items():
        moments[name]=sum(mp.mpf(str(co))*f0(nu-i,nu-j)/(4*mp.pi)**mp.mpf("1.5")
                         for (i,j),co in sp.Poly(poly,u,v).terms())
    M,P,Q,EE,FF,T,U,V,X,Y=(moments[k] for k in ("M","P","Q","E","F","T","U","V","X","Y"))
    actual=[M,P+Q*c3*c3,P+Q*c1*c1,
            2*T+(4*T-X)*e*e+2*U*(c1*c1+c3*c3)
            +(8*U-Y-EE/2)*e*c1*c3+(2*V-FF/2)*c1*c1*c3*c3]
    den=32*(nu-3)*(nu-2)*(4*nu-7)*(4*nu-5)
    l12=base*(nu-1)*(2*nu-3)/(8*(nu-2)*(4*nu-7)*(4*nu-5))
    angular=((nu-1)*(2*nu-3)*(2*nu-5)
        -(nu-1)*(2*nu-3)**2*(c1*c1+c3*c3-3*c1*c1*c3*c3)
        -2*(nu-1)*(nu-2)*(2*nu-3)*mp.sin(mp.mpf("1.4"))*mp.sin(mp.mpf("2.2"))*mp.cos(mp.mpf(".8"))
        +(nu-2)*(2*nu-5)*(2*nu-3)*mp.sin(mp.mpf(".7"))**2*mp.sin(mp.mpf("1.1"))**2*mp.cos(mp.mpf("1.6")))
    claimed=[3*l12,l12,l12,base/den*angular]
    maxerr=max(maxerr,max(abs(x-y) for x,y in zip(actual,claimed)))
    loop_values[mu]=actual
check("four loop-first scalar coefficients from independent Gamma sums",
      maxerr,mp.mpf("2e-16"))

# Existing final coefficients must follow from 8 c^2 16^(2 i mu) L_ij.
for mu in (".2",".73","1","2","5"):
    p1,p2,_,_=seeds.hard_coefficients(mu,25)
    z=1j*mp.mpf(mu);nu=mp.mpf(".5")-z
    c=mp.power(2,-1-2*z)*mp.gamma(mp.mpf(".5")-z)**2/mp.pi
    base=f0(nu,nu)/(4*mp.pi)**mp.mpf("1.5")
    common=base/(32*(nu-3)*(nu-2)*(4*nu-7)*(4*nu-5))
    actual=8*c*c*mp.power(16,2*z)*common
    expected=mp.power(2,1+4*z)*mp.gamma(mp.mpf(".5")-z)**4*common/mp.pi**2
    check(f"loop-first normalization and phase mu={mu}",abs(actual-expected),mp.mpf("1e-20"))
    if mp.mpf(mu) in loop_values:
        l11,l12,l21,l22=loop_values[mp.mpf(mu)]
        assembled=8*c*c*mp.power(16,2*z)*(l11*p1*p1+(l12+l21)*p1*p2+l22*p2*p2)
        old=seeds.clock_angular_coefficients(mu,25)
        t1,t3,phi=map(mp.mpf,(".7","1.1",".8"))
        x1,x3=mp.cos(t1),mp.cos(t3)
        reference=(old["A1"]+old["A2"]*(x1*x1+x3*x3)+old["A3"]*x1*x1*x3*x3
                   +old["B1"]*mp.sin(2*t1)*mp.sin(2*t3)*mp.cos(phi)
                   +old["D1"]*mp.sin(t1)**2*mp.sin(t3)**2*mp.cos(2*phi))
        check(f"loop-first assembled angular signal vs retained result mu={mu}",
              abs(assembled-reference)/abs(reference),mp.mpf("1e-12"))
print("Scope: local soft expansion, permutation/dummy-variable identities, scalar loop algebra and retained seed normalization; not a full unexpanded BD-box integral.",flush=True)
