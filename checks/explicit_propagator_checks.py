"""Check Section 1.3's exact Whittaker products for three Weyl propagators.
Compare (i) the existing Whittaker mode contraction and (ii) Hankel modes
with explicit helicity spinors. These are free-propagator checks, not a box
integral or a new independent derivation of the interacting amplitude.
"""
from pathlib import Path
from itertools import product
import mpmath as mp
import numpy as np
from aoki_majorana import EPSILON_LOWER as EPS, I, slash, greater, contraction


def explicit_greater(k, t1, t2, mu):
    p = np.linalg.norm(k)
    x1, x2 = -p*t1, -p*t2
    wm1, wp1, wm2, wp2 = [complex(mp.whitw(a, 1j*mu, -2j*x))
                           for x in (x1, x2) for a in (-.5, .5)]
    plus, minus = I+slash(k/p), I-slash(k/p)
    den = 4*np.sqrt(x1*x2)
    normal = (mu**2*wm1*wm2.conjugate()*plus
              + wp1*wp2.conjugate()*minus)/den
    anomalous = 1j*mu*(wm1*wp2.conjugate()*plus
                       - wp1*wm2.conjugate()*minus)@EPS/den
    conjugate = 1j*mu*EPS@(wp1*wm2.conjugate()*plus
                           - wm1*wp2.conjugate()*minus)/den
    return np.array([normal, anomalous, conjugate])


def explicit_lesser(k, t1, t2, mu):
    p = np.linalg.norm(k)
    x1, x2 = -p*t1, -p*t2
    wm1, wp1, wm2, wp2 = [complex(mp.whitw(a, 1j*mu, -2j*x))
                           for x in (x1, x2) for a in (-.5, .5)]
    normal = -(wp1.conjugate()*wp2*(I+slash(k/p))
               + mu**2*wm1.conjugate()*wm2*(I-slash(k/p)))/(4*np.sqrt(x1*x2))
    reversed_greater = explicit_greater(-k, t2, t1, mu)
    return np.array([normal, -reversed_greater[1].T, -reversed_greater[2].T])


def hankel_spinors(k, t1, t2, mu):
    p = np.linalg.norm(k)
    theta, phi = np.arccos(k[2]/p), np.arctan2(k[1], k[0])
    spinors = (np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)]),
               np.array([-np.exp(-1j*phi)*np.sin(theta/2), np.cos(theta/2)]))
    def modes(t, h):
        x = -p*t
        f = mp.sqrt(mp.pi*x)/2*mp.exp(mp.pi*mu/2)*mp.hankel1(.5-1j*mu, x)
        g = mp.sqrt(mp.pi*x)/2*mp.exp(-mp.pi*mu/2)*mp.hankel1(.5+1j*mu, x)
        phase = mp.exp(3j*mp.pi/4)/mp.sqrt(2)
        return complex(phase*(f-h*g)), complex(phase*(f+h*g))
    result = np.zeros((3, 2, 2), complex)
    for h, spin in zip((1, -1), spinors):
        u1, v1 = modes(t1, h)
        u2, v2 = modes(t2, h)
        result[0] += u1*u2.conjugate()*np.outer(spin, spin.conjugate())
        result[1] += u1*v2.conjugate()*np.outer(spin, EPS@spin.conjugate())
        result[2] += v1*u2.conjugate()*np.outer(EPS@spin, spin.conjugate())
    return result


def select(full):
    return np.array([full[:2, 2:], full[:2, :2], full[2:, 2:]])


def residual(a, b):
    return np.linalg.norm(a-b)/max(1., np.linalg.norm(b))


def main():
    lines = []
    vectors = (np.array([.3, -.4, .5]), np.array([-.8, .2, .4]), np.array([0., 0., 1.2]))
    times = ((-.31, -2.4), (-1.7, -.22), (-.8, -.8))
    previous = None
    for digits in (30, 60):
        mp.mp.dps = digits
        samples = []
        for mu in (.2, .73, 2., 5.):
            errors = dict(whittaker=0., hankel=0., sk=0., conjugation=0., canonical=0.)
            for k, (t1, t2) in product(vectors, times):
                gt, lt = explicit_greater(k, t1, t2, mu), explicit_lesser(k, t1, t2, mu)
                samples.extend((gt, lt))
                errors['whittaker'] = max(errors['whittaker'], residual(gt, select(greater(k, t1, t2, mu))))
                errors['hankel'] = max(errors['hankel'], residual(gt, hankel_spinors(k, t1, t2, mu)))
                errors['conjugation'] = max(errors['conjugation'], residual(gt[2], explicit_greater(k, t2, t1, mu)[1].conj().T))
                theta = float(t1 > t2)+.5*float(t1 == t2)
                for a, b in product((1, -1), repeat=2):
                    value = gt if (a, b) == (-1, 1) else lt if (a, b) == (1, -1) else (
                        theta*gt+(1-theta)*lt if a == 1 else (1-theta)*gt+theta*lt)
                    errors['sk'] = max(errors['sk'], residual(value, select(contraction(k, t1, t2, a, b, mu))))
                if t1 == t2:
                    errors['canonical'] = max(errors['canonical'], residual(gt[0]-lt[0], I),
                                               np.linalg.norm(gt[1]-lt[1]), np.linalg.norm(gt[2]-lt[2]))
            assert max(errors.values()) < 2e-13, (mu, errors)
            line = f'PASS mu={mu}, dps={digits}, 3 directions x 3 time pairs: '+', '.join(f'{name}={val:.3e}' for name, val in errors.items())
            print(line, flush=True); lines.append(line)
        samples = np.array(samples)
        if previous is not None:
            error = residual(samples, previous)
            assert error < 2e-13
            line = f'PASS 30 -> 60 digit stability (matrices stored in complex128): {error:.3e}'
            print(line, flush=True); lines.append(line)
        previous = samples
    lines.append('Scope: exact free propagators, all three types and all SK assignments; no late-time approximation. Not an independent box-amplitude verification.')
    Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
