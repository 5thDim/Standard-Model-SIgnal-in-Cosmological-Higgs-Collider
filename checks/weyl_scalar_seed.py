"""Top-quark hard blocks from Qin--Xianyu's branch-resolved scalar seed.

The scalar seed is analytically continued to nu = mu - i/2. Equation (92)
of arXiv:2301.07047 uses REGULARIZED 3F2, defined in its Eq. (138).
The former implementation omitted its two denominator Gamma factors.
The nu -> -nu symmetry makes all hypergeometric series convergent without
an auxiliary mass-index regulator. hard_coefficients cancels the endpoint
poles analytically before evaluation. regulated_hard_coefficients is an
independent check of the common-regulator limit, not the main evaluator.
"""

from functools import lru_cache
import mpmath as mp


def scalar_seed_folded(a, b, p, q, nu):
    """Published I_ab^{p,q}(1,1), with no factor of H**2."""
    # The scalar kernel is even in nu. Choose its sign so the regularized
    # 3F2 in the equal-SK seed has positive convergence balance.
    if a == b and mp.re(mp.mpf('.5') + a*1j*nu) < mp.mpf('.5'):
        nu = -nu
    sigma = p + q
    z = 1j * nu
    common = (mp.gamma(mp.mpf('2.5') + p - z)
              * mp.gamma(mp.mpf('2.5') + p + z)
              * mp.gamma(mp.mpf('2.5') + q - z)
              * mp.gamma(mp.mpf('2.5') + q + z)
              / (mp.gamma(3 + p) * mp.gamma(3 + q)))
    power = mp.power(2, -5 - sigma)
    if a != b:
        return mp.exp(-1j * a * mp.pi * (p - q) / 2) * power * common

    phase = mp.exp(-1j * a * mp.pi * sigma / 2)
    az = a * z
    hyp = mp.hyper([5 + sigma, mp.mpf('.5') + az, 1],
                   [mp.mpf('3.5') + p + az,
                    mp.mpf('3.5') + q + az], 1)
    particular = (mp.gamma(5 + sigma)
                  * mp.gamma(mp.mpf('2.5') + p + az)
                  * mp.gamma(mp.mpf('2.5') + q + az) * hyp
                  / (mp.gamma(mp.mpf('3.5') + p + az)
                     * mp.gamma(mp.mpf('3.5') + q + az)))
    return phase * power * (a * 1j * mp.exp(-mp.pi * nu) * common
                            - particular)


def seed_radial_derivative(a, b, p, q, nu, value):
    """d/dr I_ab(u(r),u(r)) at r=1, from the bootstrap equation."""
    p,q,nu = map(mp.mpc,(p,q,nu))
    sigma = p + q
    source = 0
    if a == b:
        source = (mp.exp(-1j * a * mp.pi * sigma / 2)
                  * mp.gamma(5 + sigma) * mp.power(2, -5 - sigma))
    bp = (nu**2 + (p + mp.mpf('2.5'))**2) / (p + 3)
    bq = (nu**2 + (q + mp.mpf('2.5'))**2) / (q + 3)
    return (bp + bq) * value / 2 - source * (1/(p + 3) + 1/(q + 3)) / 2


def regulated_hard_coefficients(mu, delta, eta, dps=40):
    """Return P1, P2, H21, H11 at a finite common endpoint regulator.

    Set eta=0 for the physical scalar index. Nonzero eta is useful only
    for checking continuation paths as both regulators vanish. At eta != 0
    the original-mu fermion IBP is only a diagnostic continuation.
    """
    with mp.workdps(dps):
        mu, delta, eta = map(mp.mpf, (mu, delta, eta))
        nu = mu - 1j/2 + 1j*eta
        base = -2 + 1j*mu + delta

        @lru_cache(None)
        def seed(a, b, p, q):
            return scalar_seed_folded(a, b, p, q, nu)

        h21 = 0
        p2 = 0
        for a in (1, -1):
            for b in (1, -1):
                for j in (0, 1):
                    for k in (0, 1):
                        p, q = base+j, base+k
                        coeff = (1j*a)**j * (1j*b)**k
                        value = seed(a, b, p, q)
                        h21 -= coeff * value
                        radial = seed_radial_derivative(a, b, p, q, nu, value)
                        p2 -= coeff * (radial - (5+p+q)*value)

        # IBP of R11 = r^{-1}(-i d_x + mu/x) R21.  With the shared
        # endpoint regulator alpha=1+i mu+delta, the coefficient
        # i(alpha-1)+mu is i*delta. Keep it until the full SK sum and
        # finite part. For eta=0 the complete H21 is finite and this term
        # ultimately vanishes, but dropping it at finite delta is incorrect.
        # The x=0 surface term cancels only after the complete SK sum.
        h11 = 0
        for a in (1, -1):
            for b in (1, -1):
                p = 1j*mu + delta
                q = -2 + 1j*mu + delta
                h11 -= 1j*(seed(a, b, p, q)
                            + 1j*b*seed(a, b, p, q+1))
        h11 += 1j*delta*h21
        p1 = h21 + 2j/(1+2j*mu)*h11
        return +p1, +p2, +h21, +h11


def regulated_diagonal_seed(mu, alpha, beta, delta, dps=40):
    """H11(alpha+delta,beta+delta;1) from the general radial IBP identity."""
    with mp.workdps(dps):
        mu, delta = mp.mpf(mu), mp.mpf(delta)
        alpha, beta = mp.mpc(alpha)+delta, mp.mpc(beta)+delta
        nu = mu - 1j/2

        @lru_cache(None)
        def seed(a, b, p, q):
            return scalar_seed_folded(a, b, p, q, nu)

        h21 = -sum((1j*a)**j*(1j*b)**k
                   * seed(a,b,alpha-3+j,beta-2+k)
                   for a in (1,-1) for b in (1,-1)
                   for j in (0,1) for k in (0,1))
        return ((1j*(alpha-1)+mu)*h21
                - 1j*sum(seed(a,b,alpha-1,beta-2)
                         + 1j*b*seed(a,b,alpha-1,beta-1)
                         for a in (1,-1) for b in (1,-1)))


def mixed_hard_seed(mu, dps=40, step='0.001'):
    """Numerical mixed H11(-i*mu,+i*mu;1); four-level Richardson limit.

    Repeat with step/2 and/or increased dps to estimate the error. This
    numerical extrapolation is separate from the exact clock coefficients.
    """
    with mp.workdps(dps):
        mu, step = mp.mpf(mu), mp.mpf(step)
        vals = [regulated_diagonal_seed(mu, -1j*mu, 1j*mu,
                                       step/2**j, dps) for j in range(4)]
        for order in range(1,4):
            vals = [(2**order*vals[j+1]-vals[j])/(2**order-1)
                    for j in range(len(vals)-1)]
        return +vals[0]


def hard_coefficients(mu, dps=40):
    """Return P1, P2, H21, H11 with the SK endpoint poles cancelled exactly.

    Only five distinct, convergent ordinary 3F2(1) values are needed.
    The balance is 1+i*mu. No mass-index regulator or extrapolation is used.
    """
    with mp.workdps(dps):
        mu = mp.mpf(mu)
        z = 1j*mu
        ch = mp.cosh(mp.pi*mu)

        @lru_cache(None)
        def q(j, k):
            return (mp.gamma(1+2*z+j+k)
                    / ((1+2*z+j)*(1+2*z+k))
                    * mp.hyper([1+2*z+j+k, 1+z, 1],
                               [2+2*z+j, 2+2*z+k], 1))

        def weight(j, k):
            return mp.power(2, -1-2*z-j-k)

        def piece(j, k):
            ans = 2*ch*weight(j, k)*q(min(j,k), max(j,k))
            if j == k == 0:
                ans -= (mp.pi**2*weight(0, 0)
                        * (mp.gamma(1+2*z)/mp.gamma(1+z))**2)
            return ans

        h21 = sum(piece(j,k) for j in (0,1) for k in (0,1))
        h11 = -2j*ch*sum(weight(2,k)*q(k,2) for k in (0,1))
        p1 = h21 + 2j/(1+2*z)*h11
        p2 = 0
        for j in (0,1):
            for k in (0,1):
                factor = ((j*(1+2*z+j)/(1+z+j)
                           + k*(1+2*z+k)/(1+z+k))/2
                          - (1+2*z+j+k))
                p2 += factor*piece(j,k)
                p2 += (ch*weight(j,k)*mp.gamma(1+2*z+j+k)
                       * (1/(1+z+j)+1/(1+z+k)))
        return +p1, +p2, +h21, +h11


def clock_angular_coefficients(mu, dps=40):
    """A1,A2,A3,B1,D1 in the note's (ks^2/(4*k12*k34))^(2i*mu) convention.

    Includes the four directed s-cut cycles and the Pauli trace. The
    overall -Nc*gt^4*H^4*ks^5/(16*prod(ki^3)*kL*kR) and c.c. are external.
    """
    with mp.workdps(dps):
        mu = mp.mpf(mu)
        v = mp.mpf('.5')-1j*mu
        p1,p2,_,_ = hard_coefficients(mu,dps)
        f0 = (mp.gamma(2*v-mp.mpf('1.5'))*mp.gamma(mp.mpf('1.5')-v)**2
              / (mp.gamma(v)**2*mp.gamma(3-2*v)))
        bubble = (mp.power(2,1+4j*mu)*mp.gamma(v)**4*f0
                  / (32*mp.pi**2*(4*mp.pi)**mp.mpf('1.5')
                     * (v-3)*(v-2)*(4*v-7)*(4*v-5)))
        a1 = bubble*(v-1)*(2*v-3)*(4*(v-3)*(3*p1**2+2*p1*p2)
                                         +(2*v-5)*p2**2)
        a2 = -bubble*p2**2*(v-1)*(2*v-3)**2
        b1 = -2*bubble*p2**2*(v-1)*(v-2)*(2*v-3)
        d1 = bubble*p2**2*(v-2)*(2*v-5)*(2*v-3)
        return dict(A1=+a1,A2=+a2,A3=-3*a2,B1=+b1,D1=+d1)


def self_check():
    """Check published branch formulae in a nonsingular real-index case."""
    with mp.workdps(30):
        nu = mp.mpf('.7')
        vals = {(a,b): scalar_seed_folded(a,b,0,0,nu)
                for a in (1,-1) for b in (1,-1)}
        assert abs(vals[1,-1]-vals[-1,1]) < mp.mpf('1e-25')
        assert abs(vals[1,1]-mp.conj(vals[-1,-1])) < mp.mpf('1e-25')
    print('PASS real-index branch symmetry')
    # An independent, elementary time integral for a conformal scalar.
    # At nu=-i/2, S^> = xy exp(i(x-y))/2 at r=1.
    with mp.workdps(30):
        for a,b in ((1,1),(-1,-1),(1,-1),(-1,1)):
            p, q = mp.mpf('.2'), mp.mpf('.7')
            if a == b:
                exact = (-mp.gamma(p+q+4)/(2*(2j*a)**(p+q+4))
                         * (1/(p+2)+1/(q+2)))
            else:
                exact = (mp.gamma(p+2)*mp.gamma(q+2)
                         / (2*(2j*a)**(p+2)*(2j*b)**(q+2)))
            assert abs(scalar_seed_folded(a,b,p,q,-.5j)-exact) < mp.mpf('1e-25')
    print('PASS conformal scalar seed against elementary time integrals')
    with mp.workdps(30):
        for a,b,expected in ((1,1,-mp.mpf(5)/8),(1,-1,mp.mpf(1)/16)):
            val=scalar_seed_folded(a,b,0,0,-.5j)
            derivative=seed_radial_derivative(a,b,0,0,-.5j,val)
            assert abs(derivative-expected) < mp.mpf('1e-25')
    print('PASS folded derivative against elementary conformal-scalar derivative')


if __name__ == '__main__':
    self_check()
    for mu in ('.2','.73','1','2','5'):
        p1,p2,_,_ = hard_coefficients(mu,30)
        print('mu=',mu,'P1=',mp.nstr(p1,16),'P2=',mp.nstr(p2,16),flush=True)
