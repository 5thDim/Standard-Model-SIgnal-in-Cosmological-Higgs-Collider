"""Independent normalization and convergence audit for the physical hard seeds.
Run from the project root: python3 checks/weyl_seed_checks.py
Writes weyl_seed_checks.txt. This does not integrate the unexpanded box.
"""
from pathlib import Path
import mpmath as mp
from weyl_scalar_seed import (
    scalar_seed_folded, hard_coefficients, regulated_hard_coefficients,
    mixed_hard_seed, clock_angular_coefficients, self_check,
)

mp.mp.dps = 35
lines = []


def check(name, residual, tol):
    assert residual < tol, (name, residual, tol)
    message = f'PASS {name}: residual {float(residual):.3e}'
    lines.append(message)
    print(message, flush=True)


self_check()
lines.append('PASS elementary conformal-scalar normalization, folded derivative, and real-index conjugation')
# This parity follows from Hankel identities, independently of the folded formula.
nu = mp.mpc('.73', '-.5')
x, y = mp.mpf('.7'), mp.mpf('1.3')
def kernel(n):
    return (mp.pi/4*mp.exp(-mp.pi*n)*(x*y)**mp.mpf('1.5')
            * mp.hankel1(1j*n,x)*mp.hankel2(-1j*n,y))
check('complex scalar-kernel nu parity', abs(kernel(nu)-kernel(-nu)), mp.mpf('1e-28'))
# With real nu neither expression is automatically sign-flipped in the code.
for a in (1,-1):
    v1 = scalar_seed_folded(a,a,mp.mpf('.2'),mp.mpf('.7'),mp.mpf('.73'))
    v2 = scalar_seed_folded(a,a,mp.mpf('.2'),mp.mpf('.7'),mp.mpf('-.73'))
    check(f'equal-SK seed nu parity, a={a}', abs(v1-v2), mp.mpf('1e-28'))

for mu in ('.2','.73','2','5'):
    low = hard_coefficients(mu,22)
    high = hard_coefficients(mu,35)
    check(f'physical hard coefficients precision, mu={mu}',
          max(abs(a-b) for a,b in zip(low,high)), mp.mpf('1e-17'))
    lines.append(f'mu={mu} P1={mp.nstr(high[0],18)} P2={mp.nstr(high[1],18)}')
    if mu == '.73':
        target = high

# Validate the analytically cancelled result against the original branch sum.
# Two-point extrapolation has O(delta^2) error; vary the continuation path too.
for ratio in (0,1,2):
    coarse = regulated_hard_coefficients('.73','0.0001',mp.mpf('0.0001')*ratio,30)
    fine = regulated_hard_coefficients('.73','0.00005',mp.mpf('0.00005')*ratio,30)
    check(f'full SK common-regulator limit, eta/delta={ratio}',
          max(abs(2*b-a-c) for a,b,c in zip(coarse,fine,target)), mp.mpf('2e-7'))

# Independently reconstruct the phase/normalization from the soft mode c,
# both Bose sums (4), the Pauli trace (2), and the scalar bubble moments.
mu=mp.mpf('.73'); z=1j*mu; v=mp.mpf('.5')-z
p1,p2=target[:2]
f0=(mp.gamma(2*v-mp.mpf('1.5'))*mp.gamma(mp.mpf('1.5')-v)**2
    /(mp.gamma(v)**2*mp.gamma(3-2*v)))
M=3*(v-1)*(2*v-3)/(8*(v-2)*(4*v-7)*(4*v-5))
P=M/3
T=(2*v-5)*(2*v-3)/(64*(v-3)*(v-2)*(4*v-7)*(4*v-5))
X=-(2*v-5)*(2*v-3)/(16*(v-2)*(4*v-7)*(4*v-5))
# theta1=theta3=pi/2, phi=pi/4 isolates A1 after azimuthal averaging.
polynomial=p1**2*M+2*p1*p2*P+p2**2*(4*T-X/2)
c=mp.power(2,-1-2*z)*mp.gamma(v)**2/mp.pi
ks,kL,kR=mp.mpf('.01'),mp.mpf('.8'),mp.mpf('1.7')
raw=8*c**2*f0/(4*mp.pi)**mp.mpf('1.5')*polynomial*(ks**2/(kL*kR))**(2*z)
coeffs=clock_angular_coefficients(mu,35)
final=coeffs['A1']*(ks**2/(16*kL*kR))**(2*z)
check('complete clock normalization and ratio phase',abs(raw-final),mp.mpf('1e-28'))
lines.append('mu=.73 angular coefficients: '+', '.join(f'{k}={mp.nstr(v,18)}' for k,v in coeffs.items()))

mixed=mixed_hard_seed('.73',30,'.001')
refined=mixed_hard_seed('.73',35,'.0005')
check('mixed seed regulator step and precision',abs(mixed-refined),mp.mpf('1e-12'))
check('mixed seed reality',abs(refined.imag),mp.mpf('1e-20'))
lines.append('mu=.73 mixed H11='+mp.nstr(refined,20))
lines.append('Scope: factorized leading nonanalytic terms; no full unexpanded-box quadrature or analytic-contact renormalization.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[-3:]))
