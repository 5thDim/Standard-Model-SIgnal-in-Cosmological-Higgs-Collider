"""Aoki two-component Majorana rules; independent Whittaker implementation.
All fields rescaled, Fourier exp(+ik.x), no i in contractions.
D denotes the paper's rescaled propagator tilde D, with type order (L,Lbar).
D^{rs}=<T_C zeta_r zeta_s^T>, zeta=(tilde psi,tilde psi^dagger_lower).
J=diag(epsilon,-epsilon), O=zeta^T J zeta/2.
"""
from functools import lru_cache
import mpmath as mp
import numpy as np

SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)
# epsilon_{12}=+1; epsilon^{alpha beta} is -EPSILON_LOWER.
EPSILON_LOWER = np.array([[0, 1], [-1, 0]], complex)
I = np.eye(2, dtype=complex)
Z = np.zeros((2, 2), complex)
J = np.block([[EPSILON_LOWER, Z], [Z, -EPSILON_LOWER]])
T = np.block([[I, Z], [Z, EPSILON_LOWER]])  # zeta=T(psi,psi^dagger_raised)
X = np.block([[Z, I], [I, Z]])
UB = np.block([[I, -I], [I, I]])/np.sqrt(2)
BETA = np.block([[I, Z], [Z, -I]])


def slash(k):
    return np.einsum('i,ijk->jk', k, SIGMA)


def uv(mu, x, helicity):
    return _uv(mu, x, helicity, mp.mp.prec)


@lru_cache(maxsize=4096)
def _uv(mu, x, helicity, precision):
    """Eq. (2.18), with mu=m/H; x=-k*tau>0."""
    low = mu/mp.sqrt(2*x)*mp.whitw(-mp.mpf('.5'), 1j*mu, -2j*x)
    high = 1j/mp.sqrt(2*x)*mp.whitw(mp.mpf('.5'), 1j*mu, -2j*x)
    return (low, high) if helicity == 1 else (high, low)


def greater(k, t1, t2, mu):
    """All four lower-index two-component Wightman contractions."""
    k = np.asarray(k, float)
    p = np.linalg.norm(k)
    blocks = [Z.copy() for _ in range(4)]
    for h in (1, -1):
        u1, v1 = map(complex, uv(mu, -p*t1, h))
        u2, v2 = map(complex, uv(mu, -p*t2, h))
        P = (I+h*slash(k/p))/2
        blocks[0] -= u1*v2.conjugate()*P@EPSILON_LOWER
        blocks[1] += u1*u2.conjugate()*P
        blocks[2] -= v1*v2.conjugate()*EPSILON_LOWER@P@EPSILON_LOWER
        blocks[3] += v1*u2.conjugate()*EPSILON_LOWER@P
    return np.block([[blocks[0], blocks[1]], [blocks[2], blocks[3]]])


def contraction(k, t1, t2, a, b, mu):
    gt = greater(k, t1, t2, mu)
    lt = -greater(-np.asarray(k), t2, t1, mu).T
    if (a, b) == (-1, 1):
        return gt
    if (a, b) == (1, -1):
        return lt
    theta = float(t1 > t2) + .5*float(t1 == t2)
    return theta*gt+(1-theta)*lt if a == 1 else (1-theta)*gt+theta*lt


def dirac_hankel(k, t1, t2, a, b, mu):
    """Independent old Dirac mode sum in mass diagonal basis."""
    p = np.linalg.norm(k)
    def fg(t):
        x = -p*t
        return [complex(mp.sqrt(mp.pi*x)/2*mp.exp(s*mp.pi*mu/2)
                        *mp.hankel1(.5-s*1j*mu, x)) for s in (1, -1)]
    f1, g1 = fg(t1)
    f2, g2 = fg(t2)
    S = slash(np.asarray(k)/p)
    U1 = np.vstack((f1*I, g1*S)); U2 = np.vstack((f2*I, g2*S))
    V1 = np.vstack((g1.conjugate()*I, -f1.conjugate()*S))
    V2 = np.vstack((g2.conjugate()*I, -f2.conjugate()*S))
    gt, lt = U1@U2.conj().T, -V1@V2.conj().T
    if (a, b) == (-1, 1): return gt
    if (a, b) == (1, -1): return lt
    theta = float(t1 > t2) + .5*float(t1 == t2)
    return theta*gt+(1-theta)*lt if a == 1 else (1-theta)*gt+theta*lt


def pairings(items):
    """Wick pairings with permutation sign, independent of loop traces."""
    if not items:
        yield 1, []
        return
    for j in range(1, len(items)):
        for sign, pairs in pairings(items[1:j]+items[j+1:]):
            yield (-1)**(j-1)*sign, [(items[0], items[j])]+pairs


def connected(pairs):
    adj = [set() for _ in range(4)]
    for i, j in pairs:
        if i//2 == j//2: return False
        adj[i//2].add(j//2); adj[j//2].add(i//2)
    reached, pending = set(), [0]
    while pending:
        v = pending.pop()
        if v not in reached:
            reached.add(v); pending.extend(adj[v]-reached)
    return len(reached) == 4


CONNECTED = [(s, p) for s, p in pairings(list(range(8))) if connected(p)]


def wick_four(D):
    """Connected <O1 O2 O3 O4> by all slots/types and 8-field Wick theorem.
    D shape (4,4,4,4): vertex,vertex,spin+type,spin+type.
    Includes each O=1/2 zeta^T J zeta; no closed-loop rule assumed.
    """
    import itertools
    slots = [(a, b, J[a, b]) for a, b in zip(*np.nonzero(J))]
    answer = 0j
    for vertices in itertools.product(slots, repeat=4):
        idx = [a for v in vertices for a in v[:2]]
        weight = np.prod([v[2]/2 for v in vertices])
        for sign, pairs in CONNECTED:
            product = sign*weight
            for i, j in pairs:
                product *= D[i//2, j//2, idx[i], idx[j]]
            answer += product
    return answer


def cycle_sum(D):
    import itertools
    return -.5*sum(np.trace(J@D[i,j]@J@D[j,k]@J@D[k,l]@J@D[l,i])
                   for j,k,l in itertools.permutations((1,2,3)) for i in (0,))
