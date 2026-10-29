"""Dirac box Wick counting, physical routing and double-soft factorization.

Run: PYTHONDONTWRITEBYTECODE=1 python3 checks/dirac_box_cut_checks.py
Direct Wick enumeration uses arbitrary 4 x 4 contractions. Physical checks
use the Whittaker modes, all 16 SK assignments and exact nonlocal branches.
No loop or time integral is evaluated in this step.
"""

from functools import lru_cache
from itertools import permutations, product
from pathlib import Path

import mpmath as mp
import numpy as np

from dirac_propagator_checks import IDENTITY, G0, spatial_gamma, wightman, select_sk


def pairings(slots):
    if not slots:
        yield 1, []
        return
    for j in range(1, len(slots)):
        for sign, pairs in pairings(slots[1:j] + slots[j+1:]):
            yield (-1)**(j-1) * sign, [(slots[0], slots[j])] + pairs


def connected_dirac_pairings():
    """Slot order is bar(T_1),T_1,...,bar(T_4),T_4."""
    result = []
    for sign, pairs in pairings(list(range(8))):
        if any(first % 2 == second % 2 for first, second in pairs):
            continue
        mapping = {}
        for first, second in pairs:
            particle, adjoint = (first, second) if first % 2 else (second, first)
            mapping[particle//2] = adjoint//2
        cycle = [0]
        while mapping[cycle[-1]] not in cycle:
            cycle.append(mapping[cycle[-1]])
        if len(cycle) == 4:
            result.append((sign, pairs, tuple(cycle)))
    return result


def wick_sum(contractions, matchings):
    result = 0j
    for indices in product(range(4), repeat=4):
        for sign, pairs, _ in matchings:
            term = complex(sign)
            for first, second in pairs:
                i, j = first//2, second//2
                if first % 2:
                    term *= contractions[i, j, indices[i], indices[j]]
                else:
                    term *= -contractions[j, i, indices[j], indices[i]]
            result += term
    return result


def endpoint(k, time, mu, eta):
    p = np.linalg.norm(k)
    f = mp.hyp0f1(.5 + 1j*eta*mu, -p*p*time*time/4)
    g = (1j*eta*p*time/(1+2j*eta*mu)
         * mp.hyp0f1(1.5+1j*eta*mu, -p*p*time*time/4))
    return complex(f)*IDENTITY + complex(g)*spatial_gamma(k/p)


def central(k, eta):
    return (IDENTITY + eta*G0) @ spatial_gamma(k/np.linalg.norm(k))


def soft_prefactor(k, time1, time2, mu, eta):
    p = np.linalg.norm(k)
    return -complex(mp.gamma(.5-1j*eta*mu)**2/(4*mp.pi)
                    * mp.power(p*p*time1*time2/4, 1j*eta*mu))


@lru_cache(maxsize=4096)
def full_line(k, time1, time2, mu, a, b):
    greater, lesser = wightman(np.array(k), time1, time2, mu)
    return select_sk(greater, lesser, time1, time2, a, b)


@lru_cache(maxsize=4096)
def cut_line(k, time1, time2, mu, eta):
    return wightman(np.array(k), time1, time2, mu, eta)[0]


def main():
    mp.mp.dps = 40
    lines = []
    def record(text):
        print(text, flush=True)
        lines.append(text)

    matchings = connected_dirac_pairings()
    cycles = {(0,)+tail for tail in permutations((1, 2, 3))}
    assert len(matchings) == 6 and {m[2] for m in matchings} == cycles
    for sign, pairs, _ in matchings:
        assert sign * (-1)**sum(first % 2 == 0 for first, _ in pairs) == -1
    rng = np.random.default_rng(20261029)
    errors = []
    for _ in range(3):
        contractions = (rng.normal(size=(4, 4, 4, 4))
                        + 1j*rng.normal(size=(4, 4, 4, 4))) / 3
        direct = wick_sum(contractions, matchings)
        traces = -sum(np.trace(contractions[a, b] @ contractions[b, c]
                               @ contractions[c, d] @ contractions[d, a])
                      for a, b, c, d in cycles)
        errors.append(abs(direct-traces)/max(1, abs(traces)))
    assert max(errors) < 2e-13
    colors = sum(all(labels[i] == labels[(i+1) % 4] for i in range(4))
                 for labels in product(range(3), repeat=4))
    assert colors == 3
    assert len(list(permutations(range(4)))) == 24
    record('PASS independent eight-field Wick enumeration: six connected directed cycles, '
           f'each with sign -1; explicit index-sum residual {max(errors):.3e}; color trace 3.')

    selected = set()
    for cycle in cycles:
        edges = {frozenset((cycle[i], cycle[(i+1) % 4])) for i in range(4)}
        if frozenset((0, 1)) in edges and frozenset((2, 3)) in edges:
            selected.add(cycle)
    rows = [(0, 1, 2, 3), (0, 1, 3, 2), (1, 0, 2, 3), (1, 0, 3, 2)]
    canonical = {row[row.index(0):]+row[:row.index(0)] for row in rows}
    assert canonical == selected

    ks = np.array([0., 0., .19])
    k_left, k_right = np.array([.71, .23, .94]), np.array([-.53, .62, .38])
    external = np.array([k_left+ks/2, -k_left+ks/2, k_right-ks/2, -k_right-ks/2])
    q = np.array([.047, -.031, .026])
    big_k, ell = q+ks, q+ks/2
    for a, b, c, d in rows:
        momenta = [q+external[a], big_k, q-external[d], q]
        for j, vertex in enumerate((a, b, c, d)):
            assert np.linalg.norm(momenta[j]-momenta[j-1]-external[vertex]) < 5e-16
        assert np.linalg.norm(momenta[0] - ((k_left if a == 0 else -k_left)+ell)) < 5e-16
        assert np.linalg.norm(momenta[2] - ((k_right if c == 2 else -k_right)+ell)) < 5e-16
    record('PASS exactly four s-channel cuts; all physical endpoint assignments and centered momentum routes conserve momentum.')

    endpoint_error = 0.
    for mu, eta, k, (t1, t2) in product(
            (.2, .73, 2., 5.), (1, -1), (q, big_k, k_left), ((-.31, -1.4), (-3.2, -2.6))):
        direct = cut_line(tuple(k), t1, t2, mu, eta)
        factored = (soft_prefactor(k, t1, t2, mu, eta)
                    * endpoint(k, t1, mu, eta) @ central(k, eta) @ endpoint(k, t2, mu, eta))
        endpoint_error = max(endpoint_error, np.linalg.norm(direct-factored)/np.linalg.norm(direct))
    assert endpoint_error < 3e-13
    record(f'PASS exact four-dimensional endpoint factorization against Whittaker branches: {endpoint_error:.3e}.')

    times = (-.43, -1.27, -.61, -1.49)
    energies = np.linalg.norm(external, axis=1)
    assert energies[0] != energies[1] and energies[2] != energies[3]
    trace_error = collected_error = weighted_error = phase_error = 0.
    for mu, eta, sk in product((.2, .73, 2.), (1, -1), product((1, -1), repeat=4)):
        prefactor = complex(mp.power(4, -2j*eta*mu)*mp.gamma(.5-1j*eta*mu)**4/(4*mp.pi)**2
                            * mp.power(np.linalg.norm(q)*np.linalg.norm(big_k), 2j*eta*mu)
                            * mp.fprod(mp.power(-t, 1j*eta*mu) for t in times))
        raw_total = 0j
        hard_blocks = {}
        scale = 0.
        for a, b, c, d in rows:
            hard_left = full_line(tuple(q+external[a]), times[a], times[b], mu, sk[a], sk[b])
            hard_right = full_line(tuple(q-external[d]), times[c], times[d], mu, sk[c], sk[d])
            soft_k = cut_line(tuple(big_k), times[b], times[c], mu, eta)
            soft_q = cut_line(tuple(q), times[d], times[a], mu, eta)
            direct = np.trace(hard_left @ soft_k @ hard_right @ soft_q)
            left = endpoint(q, times[a], mu, eta) @ hard_left @ endpoint(big_k, times[b], mu, eta)
            right = endpoint(big_k, times[c], mu, eta) @ hard_right @ endpoint(q, times[d], mu, eta)
            hard_blocks[a, b] = left
            hard_blocks[c, d] = right
            factored = prefactor*np.trace(left @ central(big_k, eta) @ right @ central(q, eta))
            norm = np.prod([np.linalg.norm(m) for m in (hard_left, soft_k, hard_right, soft_q)])
            trace_error = max(trace_error, abs(direct-factored)/norm)
            actual_phase = (soft_prefactor(big_k, times[b], times[c], mu, eta)
                            * soft_prefactor(q, times[d], times[a], mu, eta))
            phase_error = max(phase_error, abs(actual_phase-prefactor)/abs(prefactor))
            raw_total += direct
            scale += norm
        collected = prefactor*np.trace((hard_blocks[0, 1]+hard_blocks[1, 0]) @ central(big_k, eta)
                                       @ (hard_blocks[2, 3]+hard_blocks[3, 2]) @ central(q, eta))
        collected_error = max(collected_error, abs(raw_total-collected)/scale)
        # H=1: all physical external factors, vertex scale factors and SK signs.
        weight = np.prod([a*(-1/t)*(1-1j*a*k*t)*np.exp(1j*a*k*t)/(2*k**3)
                          for a, k, t in zip(sk, energies, times)])
        weighted_error = max(weighted_error, abs(weight*(raw_total-collected))/(abs(weight)*scale))
    assert max(trace_error, collected_error, weighted_error, phase_error) < 3e-13
    record('PASS 3 masses x 2 common soft branches x all 16 SK: '
           f'individual traces={trace_error:.3e}, common phase={phase_error:.3e}, '
           f'four-trace sum={collected_error:.3e}, exact unequal external weights={weighted_error:.3e}.')
    record('Scope: Dirac Wick combinatorics, momentum routing and exact cut-trace factorization. '
           'No leading-order cancellation proof, loop integral or time integral is claimed.')
    Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
