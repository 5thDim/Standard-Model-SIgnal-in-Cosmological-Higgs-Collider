"""Check the four-component modes and propagators in Section 6.

Run: PYTHONDONTWRITEBYTECODE=1 python3 checks/dirac_propagator_checks.py
The Whittaker scalar kernel uses exactly the SMS normalization D^(kappa).
Independent checks use Hankel radial solutions and exact Whittaker branches.
Matrices are stored in complex128; special functions use 30 and 50 digits.
This checks free fields and the soft expansion, not the interacting box.
"""

from functools import lru_cache
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


IDENTITY = np.eye(4, dtype=complex)
GAMMA = np.array([
    [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]],
    [[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]],
    [[0, 0, 0, -1j], [0, 0, 1j, 0], [0, 1j, 0, 0], [-1j, 0, 0, 0]],
    [[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]],
], dtype=complex)
G0 = GAMMA[0]
G5 = 1j * G0 @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
PL, PR = (IDENTITY - G5) / 2, (IDENTITY + G5) / 2
CHARGE = 1j * GAMMA[2] @ G0


def spatial_gamma(direction):
    return np.einsum('i,ijk->jk', direction, GAMMA[1:])


def polarization(direction, helicity):
    """A normalized four-component left-chiral helicity eigenvector."""
    theta = np.arccos(np.clip(direction[2], -1, 1))
    phi = np.arctan2(direction[1], direction[0])
    if helicity == 1:
        return np.array([np.cos(theta / 2),
                         np.exp(1j * phi) * np.sin(theta / 2), 0, 0])
    return np.array([-np.exp(-1j * phi) * np.sin(theta / 2),
                     np.cos(theta / 2), 0, 0])


def whittaker_branch(kappa, mu, z, eta):
    """Exact M-function branch from DLMF 13.14.2 and 13.14.33."""
    return (mp.gamma(-2j * eta * mu)
            / mp.gamma(mp.mpf('.5') - 1j * eta * mu - kappa)
            * mp.exp(-z / 2) * z ** (mp.mpf('.5') + 1j * eta * mu)
            * mp.hyp1f1(mp.mpf('.5') + 1j * eta * mu - kappa,
                        1 + 2j * eta * mu, z))


def radial(mu, x, helicity, eta=None, hankel=False):
    mu, x = mp.mpf(mu), mp.mpf(x)
    return _radial(mu, x, helicity, eta, hankel, mp.mp.prec)


@lru_cache(maxsize=4096)
def _radial(mu, x, helicity, eta, hankel, precision):
    if hankel:
        f = mp.sqrt(mp.pi * x) / 2 * mp.exp(mp.pi * mu / 2) * mp.hankel1(.5 - 1j * mu, x)
        g = mp.sqrt(mp.pi * x) / 2 * mp.exp(-mp.pi * mu / 2) * mp.hankel1(.5 + 1j * mu, x)
        phase = mp.exp(3j * mp.pi / 4) / mp.sqrt(2)
        return phase * (f - helicity * g), phase * (f + helicity * g)
    def w(kappa):
        if eta is None:
            return mp.whitw(kappa, 1j * mu, -2j * x)
        return whittaker_branch(kappa, mu, -2j * x, eta)
    low = mu / mp.sqrt(2 * x) * w(-.5)
    high = 1j / mp.sqrt(2 * x) * w(.5)
    return (low, high) if helicity == 1 else (high, low)


def positive_mode(k, time, mu, helicity, eta=None, hankel=False):
    p = np.linalg.norm(k)
    u, v = map(complex, radial(mu, -p * time, helicity, eta, hankel))
    return (u * IDENTITY + v * G0) @ polarization(k / p, helicity)


def negative_mode(k, time, mu, helicity, eta=None, hankel=False):
    u = positive_mode(k, time, mu, helicity, eta, hankel)
    return CHARGE @ G0.T @ u.conj()


def spinor_product(first, second):
    return np.outer(first, second.conj()) @ G0


def wightman(k, time1, time2, mu, eta=None, hankel=False):
    """For eta != None, select exact nonanalytic cross-branch products."""
    other = None if eta is None else -eta
    greater = sum((spinor_product(
        positive_mode(k, time1, mu, s, eta, hankel),
        positive_mode(k, time2, mu, s, other, hankel)) for s in (1, -1)))
    lesser = -sum((spinor_product(
        negative_mode(-k, time1, mu, s, other, hankel),
        negative_mode(-k, time2, mu, s, eta, hankel)) for s in (1, -1)))
    return greater, lesser


def select_sk(greater, lesser, time1, time2, a, b):
    if a != b:
        return greater if a == -1 else lesser
    theta = float(time1 > time2) + .5 * float(time1 == time2)
    if a == -1:
        theta = 1 - theta
    return theta * greater + (1 - theta) * lesser


def scalar_kernel(p, x, y, mu, kappa, lesser):
    z1, z2 = (2j * x, -2j * y) if lesser else (-2j * x, 2j * y)
    k1, k2 = (-kappa, kappa) if lesser else (kappa, -kappa)
    return (mp.exp(1j * mp.pi * kappa) / (2 * p)
            * mp.whitw(k1, 1j * mu, z1) * mp.whitw(k2, 1j * mu, z2))


def reconstructed(k, time1, time2, mu, a, b):
    """Apply the manuscript differential operators to the scalar SK kernel."""
    p = np.linalg.norm(k)
    x, y = map(mp.mpf, (-p * time1, -p * time2))
    def kernel(z, kappa):
        gt = scalar_kernel(p, z, y, mu, kappa, False)
        lt = scalar_kernel(p, z, y, mu, kappa, True)
        return select_sk(gt, lt, time1, time2, a, b)
    dm, dp = [kernel(x, kap) for kap in (-.5, .5)]
    dxm, dxp = [mp.diff(lambda z: kernel(z, kap), x) for kap in (-.5, .5)]
    rm = -x * dxm - 1j * x * dm + dm / 2
    ap = x * dxp - 1j * x * dp - dp / 2
    q = G0 @ spatial_gamma(k / p)
    plus, minus = (IDENTITY + q) / 2, (IDENTITY - q) / 2
    return float(p / mp.sqrt(x * y)) * (
        mu * (complex(dm) * minus + complex(dp) * plus)
        + 1j * (complex(rm) * plus - complex(ap) * minus) @ G0)


def soft_first_order(k, time1, time2, mu, eta):
    p = np.linalg.norm(k)
    coefficient = complex(mp.gamma(.5 - 1j * eta * mu) ** 2 / (4 * mp.pi)
                          * mp.power(p * p * time1 * time2 / 4, 1j * eta * mu))
    return -coefficient * (
        (IDENTITY + eta * G0) @ spatial_gamma(k / p)
        + 1j * p / (1 + 2j * eta * mu)
        * ((time1 - time2) * G0 - eta * (time1 + time2) * IDENTITY))


def residual(first, second):
    return np.linalg.norm(first - second) / max(1., np.linalg.norm(second))


def check_exact_algebra():
    matrices = [sp.Matrix(g).applyfunc(sp.nsimplify) for g in GAMMA]
    eye = sp.eye(4)
    metric = sp.diag(-1, 1, 1, 1)
    c = sp.I * matrices[2] * matrices[0]
    five = sp.I * matrices[0] * matrices[1] * matrices[2] * matrices[3]
    left, right = (eye - five) / 2, (eye + five) / 2
    for a, b in product(range(4), repeat=2):
        assert matrices[a] * matrices[b] + matrices[b] * matrices[a] == -2 * metric[a, b] * eye
    for g in matrices:
        assert c * g.T * c.inv() == -g
    assert c.T == -c and c.H * c == eye
    for g in matrices[1:]:
        spin = five * matrices[0] * g
        assert left * (eye + spin) / 2 + right * (eye - spin) / 2 == (eye - matrices[0] * g) / 2


def main():
    lines = []
    def record(line):
        print(line, flush=True)
        lines.append(line)

    check_exact_algebra()
    record('PASS exact Clifford, charge-conjugation and projector algebra.')
    vectors = (np.array([.3, -.4, .7]), np.array([-.8, .2, .4]), np.array([0., 0., 1.2]))
    times = ((-.31, -1.4), (-1.7, -.22), (-.8, -.8))
    masses = (.2, .73, 2., 5.)
    previous = None
    for digits in (30, 50):
        mp.mp.dps = digits
        errors = dict(hankel=0., scalar_sk=0., canonical=0., conjugation=0.,
                      charge_reversal=0., mode_equation=0., polarization=0.)
        samples = []
        for mu, k, (t1, t2) in product(masses, vectors, times):
            p = np.linalg.norm(k)
            gt, lt = wightman(k, t1, t2, mu)
            hgt, hlt = wightman(k, t1, t2, mu, hankel=True)
            samples.extend((gt, lt))
            errors['hankel'] = max(errors['hankel'], residual(gt, hgt), residual(lt, hlt))
            rgt, rlt = wightman(k, t2, t1, mu)
            errors['conjugation'] = max(errors['conjugation'],
                residual(gt.conj().T, G0 @ rgt @ G0), residual(lt.conj().T, G0 @ rlt @ G0))
            rgt, rlt = wightman(-k, t2, t1, mu)
            for a, b in product((1, -1), repeat=2):
                value = select_sk(gt, lt, t1, t2, a, b)
                reverse = select_sk(rgt, rlt, t2, t1, b, a)
                errors['scalar_sk'] = max(errors['scalar_sk'], residual(value, reconstructed(k, t1, t2, mu, a, b)))
                errors['charge_reversal'] = max(errors['charge_reversal'], residual(value, CHARGE @ reverse.T @ CHARGE.conj().T))
            if t1 == t2:
                errors['canonical'] = max(errors['canonical'], residual(gt - lt, G0))
                for kap in (-.5, .5):
                    jump = scalar_kernel(p, -p*t1, -p*t1, mu, kap, False) - scalar_kernel(p, -p*t1, -p*t1, mu, kap, True)
                    assert abs(jump) < mp.mpf('1e-25')
            for s in (1, -1):
                h = polarization(k / p, s)
                ps = (IDENTITY + s * G5 @ G0 @ spatial_gamma(k / p)) / 2
                errors['polarization'] = max(errors['polarization'], residual(np.outer(h, h.conj()), PL @ ps))
                u, v = radial(mu, -p*t1, s)
                du, dv = [mp.diff(lambda t: radial(mu, -p*t, s)[j], mp.mpf(t1)) for j in (0, 1)]
                first = 1j*du + s*p*u + mu/t1*v
                second = 1j*dv - s*p*v + mu/t1*u
                errors['mode_equation'] = max(errors['mode_equation'], float(abs(first)), float(abs(second)))
                mode_derivative = (complex(du) * IDENTITY + complex(dv) * G0) @ h
                v_derivative = CHARGE @ G0.T @ mode_derivative.conj()
                negative = negative_mode(k, t1, mu, s)
                errors['mode_equation'] = max(errors['mode_equation'], np.linalg.norm(
                    1j*v_derivative - (-p*G0@spatial_gamma(k/p) - mu/t1*G0) @ negative))
        assert max(errors.values()) < 3e-13, errors
        record(f'PASS {digits} digits, 4 masses x 3 directions x 3 time pairs, all 4 SK: '
               + ', '.join(f'{name}={value:.3e}' for name, value in errors.items()))
        samples = np.array(samples)
        if previous is not None:
            error = residual(samples, previous)
            assert error < 3e-13
            record(f'PASS 30 -> 50 digit stability (complex128 matrices): {error:.3e}')
        previous = samples

    soft_equality = 0.
    ratios = []
    for mu, k, (t1, t2), eta in product(masses, vectors[:2], times[:2], (1, -1)):
        errors = []
        for scale in (.02, .01, .005):
            gt, lt = wightman(scale*k, t1, t2, mu, eta)
            norm = np.linalg.norm(gt)
            soft_equality = max(soft_equality, np.linalg.norm(gt-lt)/norm)
            errors.append(np.linalg.norm(gt-soft_first_order(scale*k, t1, t2, mu, eta))/norm)
        ratios.extend((errors[0]/errors[1], errors[1]/errors[2]))
    assert soft_equality < 3e-13
    assert min(ratios) > 3.95 and max(ratios) < 4.05, (min(ratios), max(ratios))
    record(f'PASS exact nonlocal branches: greater=lesser residual {soft_equality:.3e}; '
           f'first-order remainder step-halving ratios {min(ratios):.6f}--{max(ratios):.6f}.')
    record('Scope: free Dirac fields, scalar reconstruction, canonical jump and soft branches; no box integral or diagram-counting validation.')
    Path(__file__).with_suffix('.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
