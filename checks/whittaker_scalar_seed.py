"""Chemical-potential scalar seeds for the two-component top box.

Direct Whittaker route, without the auxiliary Hankel modes. The scalar
propagator is normalized as Qin--Xianyu (2301.07047), Eq. (96), with
nu=mu and kappa=i*h*chemical_potential=+/-1/2. The physical fermion
has no chemical potential. Eqs. (130),(131) are continued analytically.
This is the sole scalar-seed implementation used by the top-box checks.
"""
from functools import lru_cache
import mpmath as mp


def hyper3_regularized(upper, lower):
    """3F2(1)/[Gamma(d) Gamma(e)], with a convergent Thomae representation.

    Euler transformation inside the beta integral gives a new balance e-c.
    This avoids the zero-real-balance series at kappa=+/-1/2.
    """
    upper, lower = list(upper), list(lower)
    balance = sum(lower)-sum(upper)
    choices = [(mp.re(balance), None)]
    for i in range(3):
        for j in range(2):
            choices.append((mp.re(lower[j]-upper[i]), (i,j)))
    _, choice = max(choices, key=lambda item:item[0])
    if choice is None:
        return mp.hyper(upper,lower,1)*mp.rgamma(lower[0])*mp.rgamma(lower[1])
    i,j=choice
    c=upper[i]; a,b=[v for k,v in enumerate(upper) if k!=i]
    e=lower[j]; d=lower[1-j]
    return (mp.gamma(balance)*mp.rgamma(e-c)*mp.rgamma(balance+c)
            *mp.rgamma(d)*mp.hyper([d-a,d-b,c],[d,balance+c],1))


def whittaker_seed_folded(a,b,p,q,mu,kappa):
    """Published branch-resolved I_ab^(kappa);p,q(1,1), Eqs. (130),(131)."""
    p,q,mu,kappa=map(mp.mpc,(p,q,mu,kappa))
    if a==b and (mp.re(q),mp.im(q))<(mp.re(p),mp.im(p)):
        p,q=q,p
    return _whittaker_seed_folded(a,b,p,q,mu,kappa,mp.mp.prec)


@lru_cache(maxsize=4096)
def _whittaker_seed_folded(a,b,p,q,mu,kappa,precision):
    s=p+q; z=1j*mu; t=mp.mpf('1.5')
    power=mp.power(2,-3-s)
    if a!=b:
        return (power*mp.exp(-1j*a*mp.pi*(p-q)/2+1j*mp.pi*kappa)
                *mp.gamma(t+p-z)*mp.gamma(t+p+z)
                *mp.gamma(t+q-z)*mp.gamma(t+q+z)
                *mp.rgamma(2+p+a*kappa)*mp.rgamma(2+q+b*kappa))
    phase=mp.exp(-1j*a*mp.pi*s/2)
    outside=-a*1j*phase*power*mp.gamma(t+p+a*z)*mp.gamma(t+q+a*z)
    homogeneous=((mp.exp(2j*mp.pi*kappa)+mp.exp(-2*mp.pi*mu))/(2*mp.pi)
                 *mp.gamma(t+p-a*z)*mp.gamma(t+q-a*z)
                 *mp.gamma(mp.mpf('.5')+a*kappa-z)
                 *mp.gamma(mp.mpf('.5')+a*kappa+z)
                 *mp.rgamma(2+p+a*kappa)*mp.rgamma(2+q+a*kappa))
    particular=(a*1j*mp.gamma(3+s)
                *hyper3_regularized([3+s,mp.mpf('.5')-a*kappa+a*z,1],
                                   [mp.mpf('2.5')+p+a*z,mp.mpf('2.5')+q+a*z]))
    return outside*(homogeneous+particular)


def seed_radial_derivative(a,b,p,q,mu,kappa,value):
    """d/dr I(u(r),u(r)) at r=1 from the Whittaker bootstrap equation."""
    p,q,mu,kappa=map(mp.mpc,(p,q,mu,kappa))
    source=(mp.exp(-1j*a*mp.pi*(p+q)/2)*mp.power(2,-3-p-q)
            *mp.gamma(3+p+q)) if a==b else 0
    dp,dq=2+p+a*kappa,2+q+b*kappa
    return (((mu**2+(p+mp.mpf('1.5'))**2)/dp
             +(mu**2+(q+mp.mpf('1.5'))**2)/dq)*value
            +(1/dp+1/dq)*source)/2


def hard_integrals(mu,alpha,beta,dps=40,derivative=False):
    """H0(alpha,beta;1), H1(alpha,beta;1), optionally d_r H0.

    Arguments include the common endpoint regulator. Integration by parts
    is performed first in a convergent domain, then continued as a complete
    SK sum. No endpoint terms are discarded on an individual branch.
    """
    with mp.workdps(dps):
        mu=mp.mpf(mu);alpha,beta=map(mp.mpc,(alpha,beta))
        half=mp.mpf('.5')
        @lru_cache(None)
        def moment(a,b,p,q,kappa,dr=False):
            val=whittaker_seed_folded(a,b,p,q,mu,kappa)
            if dr:
                val=(seed_radial_derivative(a,b,p,q,mu,kappa,val)
                     -(3+p+q)*val)
            return -a*b*val
        h0=h1=dh0=mp.mpc(0)
        for a in (1,-1):
            for b in (1,-1):
                for j in (0,1):
                    for l in (0,1):
                        p,q=alpha-mp.mpf('1.5')+j,beta-mp.mpf('1.5')+l
                        weight=1j/2*a*b*(1j*a)**j*(1j*b)**l
                        jp,jm=moment(a,b,p,q,half),moment(a,b,p,q,-half)
                        kp=moment(a,b,p+1,q,half) if a!=1 or derivative else 0
                        km=moment(a,b,p+1,q,-half) if a!=-1 or derivative else 0
                        h0+=weight*((p+mp.mpf('1.5')+1j*mu)*(jm-jp)
                                    +1j*(a-1)*kp-1j*(a+1)*km)
                        h1+=weight*((p+mp.mpf('1.5')-1j*mu)*(jm+jp)
                                    -1j*(a-1)*kp-1j*(a+1)*km)
                        if derivative:
                            djp,djm=(moment(a,b,p,q,half,True),
                                     moment(a,b,p,q,-half,True))
                            dkp,dkm=(moment(a,b,p+1,q,half,True),
                                     moment(a,b,p+1,q,-half,True))
                            dh0+=weight*((p+mp.mpf('1.5')+1j*mu)*(djm-djp)
                                         -1j*(kp+km)+1j*(a-1)*dkp-1j*(a+1)*dkm)
        return (+h0,+h1,+dh0) if derivative else (+h0,+h1)


def regulated_hard_coefficients(mu,delta,dps=40):
    with mp.workdps(dps):
        mu,delta=mp.mpf(mu),mp.mpf(delta);z=1j*mu
        h0,_,dh0=hard_integrals(mu,z+delta,z+delta,dps,True)
        _,h1=hard_integrals(mu,1+z+delta,z+delta,dps)
        return +(h0+2j/(1+2*z)*h1),+(dh0-h0),+h0,+h1


def _limit(evaluator,step,levels):
    values=[evaluator(step/2**j) for j in range(levels)]
    for order in range(1,levels):
        values=[tuple((2**order*y-x)/(2**order-1) for x,y in zip(a,b))
                for a,b in zip(values,values[1:])]
    return values[0]


def extrapolated_hard_coefficients(mu,dps=35,step='0.001',levels=4):
    """Numerical common-regulator limit of the Whittaker seed sums."""
    with mp.workdps(dps):
        return _limit(lambda d:regulated_hard_coefficients(mu,d,dps),
                      mp.mpf(step),levels)


def hard_coefficients(mu,dps=35):
    """Finite reduction of the complete Whittaker SK sum, no extrapolation.

    The Q representation is checked against regulated_hard_coefficients;
    it is algebraically the same physical integral, not a new seed.
    Returns P1, P2, H0, H1.
    """
    with mp.workdps(dps):
        mu=mp.mpf(mu);z=1j*mu;ch=mp.cosh(mp.pi*mu)
        @lru_cache(None)
        def Q(j,l):
            return (mp.power(2,-1-2*z-j-l)*mp.gamma(1+2*z+j+l)
                    /((1+2*z+j)*(1+2*z+l))
                    *mp.hyper([1+2*z+j+l,1+z,1],[2+2*z+j,2+2*z+l],1))
        end=mp.pi**2*mp.power(2,-1-2*z)*(mp.gamma(1+2*z)/mp.gamma(1+z))**2
        h0=2*ch*(Q(0,0)+2*Q(0,1)+Q(1,1))-end
        h1=-2j*ch*(Q(0,2)+Q(1,2))
        p2=2*ch*sum((j*(1+2*z+j)/(2*(1+z+j))
                    +l*(1+2*z+l)/(2*(1+z+l))-(1+2*z+j+l))*Q(j,l)
                   for j in (0,1) for l in (0,1))
        p2+=(1+2*z)*end
        p2+=mp.power(2,-1-2*z)*ch*mp.gamma(1+2*z)*(4+5*z+2*z*z)/(1+z)
        return +(h0+2j/(1+2*z)*h1),+p2,+h0,+h1


def mixed_hard_seed(mu,dps=35,step='0.001',levels=4):
    with mp.workdps(dps):
        mu=mp.mpf(mu)
        return _limit(lambda d:(hard_integrals(mu,-1j*mu+d,1j*mu+d,dps)[1],),
                      mp.mpf(step),levels)[0]


def clock_angular_coefficients(mu,dps=35):
    """Unchanged central bubble, using only Whittaker hard coefficients."""
    with mp.workdps(dps):
        mu=mp.mpf(mu);v=mp.mpf('.5')-1j*mu
        p1,p2,_,_=hard_coefficients(mu,dps)
        moment=(mp.gamma(2*v-mp.mpf('1.5'))*mp.gamma(mp.mpf('1.5')-v)**2
                /(mp.gamma(v)**2*mp.gamma(3-2*v)))
        bubble=(mp.power(2,1+4j*mu)*mp.gamma(v)**4*moment
                /(32*mp.pi**2*(4*mp.pi)**mp.mpf('1.5')
                  *(v-3)*(v-2)*(4*v-7)*(4*v-5)))
        a1=bubble*(v-1)*(2*v-3)*(4*(v-3)*(3*p1**2+2*p1*p2)+(2*v-5)*p2**2)
        a2=-bubble*p2**2*(v-1)*(2*v-3)**2
        b1=-2*bubble*p2**2*(v-1)*(v-2)*(2*v-3)
        d1=bubble*p2**2*(v-2)*(2*v-5)*(2*v-3)
        return dict(A1=+a1,A2=+a2,A3=-3*a2,B1=+b1,D1=+d1)
