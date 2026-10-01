"""Independent Whittaker, slot-Wick and Dirac correspondence tests.
Run from project root: python3 checks/aoki_majorana_checks.py.
No scalar-seed formula is used in these tests.
"""
from itertools import product, permutations
from pathlib import Path
import mpmath as mp
import numpy as np
from aoki_majorana import (uv, greater, contraction, dirac_hankel, EPSILON_LOWER, I, J,
                           T, X, UB, BETA, CONNECTED, wick_four, cycle_sum)
mp.mp.dps = 40
lines = []
def check(name, residual, tol):
    assert residual < tol, (name, residual, tol)
    line = f'PASS {name}: residual {float(residual):.3e}'
    print(line, flush=True); lines.append(line)

for mu in map(mp.mpf, ('.2', '.73', '2', '5')):
    residuals = []
    for x in map(mp.mpf, ('.003', '.7', '12')):
        f = mp.sqrt(mp.pi*x)/2*mp.exp(mp.pi*mu/2)*mp.hankel1(.5-1j*mu, x)
        g = mp.sqrt(mp.pi*x)/2*mp.exp(-mp.pi*mu/2)*mp.hankel1(.5+1j*mu, x)
        phase = mp.exp(3j*mp.pi/4)
        for h in (1, -1):
            u, v = uv(mu, x, h)
            residuals += [abs(u-phase*(f-h*g)/mp.sqrt(2)),
                          abs(v-phase*(f+h*g)/mp.sqrt(2)),
                          abs(abs(u)**2+abs(v)**2-1),
                          abs(-1j*mp.diff(lambda xx: uv(mu, xx, h)[0],x)+h*u-mu/x*v),
                          abs(-1j*mp.diff(lambda xx: uv(mu, xx, h)[1],x)-h*v-mu/x*u)]
    check(f'Whittaker/Hankel, first-order EOM, norm mu={mu}', max(residuals), mp.mpf('1e-32'))
    # BD condition checked by large-x convergence, separately from mode equality.
    errors = []
    for x in (mp.mpf('5000'), mp.mpf('10000')):
        u,v = uv(mu,x,1)
        errors.append(abs(v/(mp.exp(1j*mp.pi/4)*mp.exp(1j*x))-1)+abs(u))
    assert errors[1] < .51*errors[0]
    lines.append(f'PASS BD asymptotic convergence mu={mu}: ratio {float(errors[1]/errors[0]):.6f}')

rng = np.random.default_rng(20260929)
for mu in (.2, .73, 2.):
    residuals = []
    for trial in range(3):
        k = rng.normal(size=3); t1,t2 = -rng.uniform(.05,3,2)
        for a,b in product((1,-1), repeat=2):
            D = contraction(k,t1,t2,a,b,mu)
            R = UB.conj().T@T.conj().T@D@X@T@UB
            old = dirac_hankel(k,t1,t2,a,b,mu)
            residuals += [np.linalg.norm(R-old), np.linalg.norm(D+contraction(-k,t2,t1,b,a,mu).T)]
        residuals += [np.linalg.norm(greater(k,t1,t2,mu)[2:,2:]
                                      - greater(k,t2,t1,mu)[:2,:2].conj().T)]
        gt = greater(k,t1,t1,mu)
        lt = -greater(-k,t1,t1,mu).T
        residuals += [np.linalg.norm(gt-lt-X)]
    check(f'all four two-component SK blocks and canonical anticommutator mu={mu}',max(residuals),3e-13)

# Compare the full analytic 0F1 soft branch to the branch extracted from
# the Bessel decomposition of the Hankel modes, without a soft truncation.
for mu in map(mp.mpf, ('.2', '.73', '2')):
    x,y=mp.mpf('.27'),mp.mpf('.49');z=1j*mu
    def pieces(t):
        nu=.5-z;ng=mp.sqrt(mp.pi)/2*mp.exp(-mp.pi*mu/2)
        nf=mp.sqrt(mp.pi)/2*mp.exp(mp.pi*mu/2)
        fp=nf*mp.sqrt(t)*mp.besselj(-nu,t)/(1j*mp.sin(mp.pi*nu))
        fm=-nf*mp.sqrt(t)*mp.exp(-1j*mp.pi*nu)*mp.besselj(nu,t)/(1j*mp.sin(mp.pi*nu))
        nu=.5+z
        gm=ng*mp.sqrt(t)*mp.besselj(-nu,t)/(1j*mp.sin(mp.pi*nu))
        gp=-ng*mp.sqrt(t)*mp.exp(-1j*mp.pi*nu)*mp.besselj(nu,t)/(1j*mp.sin(mp.pi*nu))
        return fp,fm,gp,gm
    fp,_,gp,_=pieces(x);_,fm,_,gm=pieces(y)
    direct=mp.matrix([fp,-gp])*mp.matrix([[mp.conj(fm),mp.conj(gm)]])
    f=lambda t:mp.hyp0f1(.5+z,-t*t/4)
    g=lambda t:-1j*t/(1+2*z)*mp.hyp0f1(1.5+z,-t*t/4)
    c=mp.power(2,-1-2*z)*mp.gamma(.5-z)**2/mp.pi
    projected=c*(x*y)**z*mp.matrix([f(x),g(x)])*mp.matrix([[g(y),f(y)]])
    check(f'exact nonanalytic 0F1 soft branch vs Hankel-Bessel split mu={mu}',mp.norm(direct-projected),mp.mpf('1e-32'))

# 105 Wick matchings -> exactly 48 connected cycles: 16 slot assignments
# per undirected box. Vertex factors (1/2)^4 are retained in wick_four.
assert len(CONNECTED) == 48
lines.append('PASS Majorana combinatorics: 48 connected slot pairings = 3 undirected cycles x 16; each vertex retains 1/2')
for mu in (.2, .73, 2.):
    residuals=[]
    ts=-rng.uniform(.2,2,4)
    momenta={(i,j):rng.normal(size=3) for i in range(4) for j in range(i+1,4)}
    for branches in product((1,-1),repeat=4):
        D=np.zeros((4,4,4,4),complex); R=np.zeros_like(D)
        for (i,j),k in momenta.items():
            D[i,j]=contraction(k,ts[i],ts[j],branches[i],branches[j],mu)
            D[j,i]=-D[i,j].T
            R[i,j]=dirac_hankel(k,ts[i],ts[j],branches[i],branches[j],mu)
            R[j,i]=dirac_hankel(-k,ts[j],ts[i],branches[j],branches[i],mu)
        wick=wick_four(D)
        trace=cycle_sum(D)
        old=-sum(np.trace(BETA@R[0,j]@BETA@R[j,k]@BETA@R[k,l]@BETA@R[l,0])
                 for j,k,l in permutations((1,2,3)))
        residuals += [abs(wick-trace), abs(2*wick-old)]
    check(f'8-field Wick enumeration vs all six cycles, two Majoranas vs Dirac, 16 SK choices mu={mu}', max(residuals),2e-12)
lines.append('Scope: exact finite-time contractions on generic edge momenta; no unexpanded-box time/momentum quadrature.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
