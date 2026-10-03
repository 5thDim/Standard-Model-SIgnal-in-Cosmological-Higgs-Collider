"""Checks of the scalar-chemical-potential reconstruction.

Run: PYTHONDONTWRITEBYTECODE=1 python3 -u checks/whittaker_seed_checks.py
The old Hankel implementation is an independent comparison only.
No unexpanded-box integration is claimed.
"""
from itertools import product
import mpmath as mp
import numpy as np
from aoki_majorana import contraction, slash, EPSILON_LOWER as EPS
from aoki_weyl import scalar_coefficients
import whittaker_scalar_seed as new
import weyl_scalar_seed as old


def check(name,error,tol):
    assert error<tol,(name,error,tol)
    print(f"PASS {name}: residual {float(error):.3e}",flush=True)


def scalar(mu,kappa,X,Y,a,b):
    """Dimensionless k Dsc; differentiation is with respect to X."""
    norm=mp.exp(1j*mp.pi*kappa)/2
    gt=norm*mp.whitw(kappa,1j*mu,-2j*X)*mp.whitw(-kappa,1j*mu,2j*Y)
    lt=norm*mp.whitw(-kappa,1j*mu,2j*X)*mp.whitw(kappa,1j*mu,-2j*Y)
    if (a,b)==(-1,1): return gt
    if (a,b)==(1,-1): return lt
    if X==Y:return (gt+lt)/2
    return gt if ((X<Y)==(a==1)) else lt


mp.mp.dps=28
errors=[];coeff_errors=[];reversal_errors=[]
n=np.array([.3,-.4,.5]);n/=np.linalg.norm(n)
Pp=(np.eye(2)+slash(n))/2;Pm=(np.eye(2)-slash(n))/2
for mu in (.2,.73,2.,5.):
    for X,Y in ((.37,1.41),(1.41,.37)):
        for a,b in product((1,-1),repeat=2):
            sp,sm=(scalar(mu,k,X,Y,a,b) for k in (.5,-.5))
            ap=X*mp.diff(lambda x:scalar(mu,.5,x,Y,a,b),X)+(-1j*X-.5)*sp
            rm=-X*mp.diff(lambda x:scalar(mu,-.5,x,Y,a,b),X)+(-1j*X+.5)*sm
            root=np.sqrt(X*Y)
            normal=1j/root*(complex(rm)*Pm-complex(ap)*Pp)
            anomalous=-mu/root*(complex(sm)*Pp+complex(sp)*Pm)@EPS
            conjugate=mu/root*EPS@(complex(sp)*Pp+complex(sm)*Pm)
            exact=contraction(n,-X,-Y,a,b,mu)
            errors.extend([np.linalg.norm(normal-exact[:2,2:]),
                           np.linalg.norm(anomalous-exact[:2,:2]),
                           np.linalg.norm(conjugate-exact[2:,2:])])
            k0=1j/(2*root)*(ap-1j*mu*sp+rm+1j*mu*sm)
            k1=1j/(2*root)*(rm-1j*mu*sm-ap-1j*mu*sp)
            ref=scalar_coefficients(exact,n)
            coeff_errors.extend([abs(complex(k0)-ref[0]),abs(complex(k1)-ref[1])])
            reverse=contraction(n,-Y,-X,b,a,mu)
            rev=scalar_coefficients(reverse,n)
            reversal_errors.extend([abs(ref[0]-rev[0]),abs(ref[2]+rev[1])])
check("three Dsc-derived propagators, all SK, four masses",max(errors),2e-12)
check("two hard coefficients from direct vertex-type contraction",max(coeff_errors),2e-12)
check("hard-kernel reversal identities",max(reversal_errors),2e-12)

errors=[]
for mu in (.2,.73,2.):
    for kappa in (-.5,.5):
        for x in (.3,1.4):
            errors.append(abs(scalar(mu,kappa,x,x,-1,1)
                              -scalar(mu,kappa,x,x,1,-1)))
        remainder=[]
        for scale in (.004,.002,.001):
            z=-2j*scale
            approx=sum(mp.gamma(-2j*d*mu)/mp.gamma(.5-1j*d*mu-kappa)
                       *z**(.5+1j*d*mu)*(1-kappa*z/(1+2j*d*mu))
                       for d in (1,-1))
            remainder.append(abs(approx-mp.whitw(kappa,1j*mu,z))/mp.sqrt(abs(z)))
        assert remainder[1]<.4*remainder[0] and remainder[2]<.4*remainder[1]
check("scalar equal-time continuity",max(errors),mp.mpf('1e-23'))
print("PASS Whittaker late-time leading and subleading expansion, three masses",flush=True)

mp.mp.dps=20
errors=[];derivative_errors=[]
for a,b in product((1,-1),repeat=2):
    p,q,mu=mp.mpf('.2'),mp.mpf('.7'),mp.mpf('.73')
    val=new.whittaker_seed_folded(a,b,p,q,mu,0)
    expected=old.scalar_seed_folded(a,b,p-1,q-1,mu)
    errors.append(abs(val-expected))
    derivative_errors.append(abs(new.seed_radial_derivative(a,b,p,q,mu,0,val)
                                 -old.seed_radial_derivative(a,b,p-1,q-1,mu,expected)))
check("zero-chemical-potential seed normalization",max(errors),mp.mpf('1e-16'))
check("Whittaker bootstrap derivative at zero chemical potential",max(derivative_errors),mp.mpf('1e-16'))

for mu in ('.2','.73','2'):
    direct=new.regulated_hard_coefficients(mu,'.8',20)
    reference=old.regulated_hard_coefficients(mu,'.8',0,20)
    check(f"direct Whittaker vs independent hard seed, delta=.8, mu={mu}",
          max(abs(a-b) for a,b in zip(direct,reference)),mp.mpf('1e-14'))

# The direct Whittaker seed, not the finite-Q evaluator, is used here.
step=mp.mpf('.001')
fine=new.extrapolated_hard_coefficients('.73',25,step,4)
finite=new.hard_coefficients('.73',25)
check("Whittaker full SK regulator limit vs finite reduction, mu=.73",
      max(abs(a-b) for a,b in zip(fine,finite)),mp.mpf('2e-11'))
refined=new.extrapolated_hard_coefficients('.73',25,step/2,4)
check("clock regulator step halving",max(abs(a-b) for a,b in zip(fine,refined)),
      mp.mpf('2e-11'))
for mu in ('.2','.73','1','2','5'):
    current=new.hard_coefficients(mu,25)
    expected=old.hard_coefficients(mu,25)
    check(f"finite coefficients retained, mu={mu}",
          max(abs(a-b) for a,b in zip(current,expected)),mp.mpf('1e-20'))
    print(f"mu={mu}: P1={mp.nstr(current[0],16)}, P2={mp.nstr(current[1],16)}",flush=True)

mixed=new.mixed_hard_seed('.73',25,'.001')
refined=new.mixed_hard_seed('.73',25,'.0005')
reference=old.mixed_hard_seed('.73',25,'.0005')
check("mixed Whittaker regulator step halving",abs(mixed-refined),mp.mpf('2e-11'))
check("mixed Whittaker vs independent seed",abs(refined-reference),mp.mpf('2e-11'))
check("mixed seed reality",abs(refined.imag),mp.mpf('1e-14'))
print("mixed mu=.73: "+mp.nstr(refined,18),flush=True)
print("Scope: propagator/type identities, published seed reduction and regulator checks; "
      "not a full unexpanded BD-box quadrature.",flush=True)
