"""Soft endpoints and hard contractions in the original (psi, psi^dagger) types.
No mass-basis transformation is used. The 4x4 arrays only pack the four
lower-index 2x2 contractions; endpoints perform the explicit type sum.
"""
import mpmath as mp
import numpy as np
from aoki_majorana import I, EPSILON_LOWER as EPS, J, slash, uv


def endpoints(mu, k, tau, branch=1):
    """Return stacked L_r and concatenated U_s for either nonanalytic branch."""
    k = np.asarray(k, float)
    p = np.linalg.norm(k)
    x = -p*tau
    sig = slash(k/p)
    f = complex(mp.hyp0f1(.5+1j*mu, -x*x/4))
    g = complex(-1j*x/(1+2j*mu)*mp.hyp0f1(1.5+1j*mu, -x*x/4))
    if branch == 1:
        left = np.vstack((f*I+g*sig, EPS@(f*I-g*sig)))
        right = np.hstack(((f*I+g*sig)@EPS, f*I-g*sig))
    elif branch == -1:
        f, g = f.conjugate(), g.conjugate()
        left = np.vstack((f*I-g*sig, -EPS@(f*I+g*sig)))
        right = np.hstack(((f*I-g*sig)@EPS, -(f*I+g*sig)))
    else:
        raise ValueError('branch must be +1 or -1')
    return left, right


def soft_branch(mu, k, tau1, tau2, branch=1):
    p = np.linalg.norm(k)
    c = complex(2**(-1-2j*mu)*mp.gamma(.5-1j*mu)**2/mp.pi)
    left, _ = endpoints(mu, k, tau1, branch)
    _, right = endpoints(mu, k, tau2, branch)
    prefactor = -c/2 if branch == 1 else c.conjugate()/2
    return prefactor*(p*p*tau1*tau2)**(branch*1j*mu)*left@slash(k/p)@right


def hard_block(D, mu, q, K, tau1, tau2, branches=(1, 1)):
    """-1/2 sum_rs U_r(q) J_r D^{rs} J_s L_s(K)."""
    _, right = endpoints(mu, q, tau1, branches[0])
    left, _ = endpoints(mu, K, tau2, branches[1])
    return -.5*right@J@D@J@left


def greater_grid(mu, k, times):
    """All field contractions on tau=-times, from helicity mode sums."""
    p = np.linalg.norm(k)
    D = np.zeros((len(times), len(times), 4, 4), complex)
    for h in (1, -1):
        proj = (I+h*slash(k/p))/2
        vals = np.array([[complex(a) for a in uv(mu, float(p*t), h)]
                         for t in times])
        u, v = vals.T
        for (r, s), coef, mat in [
                ((0, 0), -np.outer(u, v.conj()), proj@EPS),
                ((0, 1), np.outer(u, u.conj()), proj),
                ((1, 0), -np.outer(v, v.conj()), EPS@proj@EPS),
                ((1, 1), np.outer(v, u.conj()), EPS@proj)]:
            D[:, :, 2*r:2*r+2, 2*s:2*s+2] += coef[:, :, None, None]*mat
    return D


def summed_grid(mu, k, e1, e2, times, weights):
    """Complete SK sum with Higgs modes and same-branch soft time powers.
    Restricted to a finite real-time interval; not a BD contour integral.
    """
    gt = greater_grid(mu, k, times)
    lt = -greater_grid(mu, -np.asarray(k), times).swapaxes(0, 1).swapaxes(2, 3)
    x, y = times[:, None], times[None, :]
    theta = ((x < y)+.5*(x == y))[:, :, None, None]
    out = np.zeros_like(gt)
    for a in (1, -1):
        for b in (1, -1):
            ex = (1+1j*a*e1*x)*np.exp(-1j*a*e1*x)*(1+1j*b*e2*y)*np.exp(-1j*b*e2*y)
            D = (gt if (a, b) == (-1, 1) else lt if (a, b) == (1, -1)
                 else theta*gt+(1-theta)*lt if a == 1 else (1-theta)*gt+theta*lt)
            out += a*b*ex[:, :, None, None]*D
    w = weights[:, None]*weights[None, :]*(x*y)**(1j*mu-1)
    return out*w[:, :, None, None]


def scalar_coefficients(D, direction):
    """GF, FF, GG coefficients after the four-type contraction, not a basis map.
    Accepts a single D or arrays with arbitrary leading dimensions.
    """
    sig = slash(direction)
    left0 = np.vstack((I, EPS))
    right0 = np.hstack((EPS, I))
    left1 = np.vstack((sig, -EPS@sig))
    right1 = np.hstack((sig@EPS, -sig))
    v0 = -.5*right0@J@D@J@left0
    vq = -.5*right1@J@D@J@left0
    vk = -.5*right0@J@D@J@left1
    trace = lambda a: np.trace(a, axis1=-2, axis2=-1)/2
    return -trace(sig@v0), trace(sig@vq), -trace(vk@sig)
