"""Counterexamples to specific printed v7_24 formulas, without fitting our box.
Sources: v7_24 (9)-(21), You v1/v2 (2.17),(2.19),(2.30).
The SK table uses our existing scalar-seed evaluator; it illustrates why
branch selection is not a regulated finite answer, rather than providing an
independent proof of that evaluator. No claim to audit the entire You paper.
"""
from pathlib import Path
from itertools import product
import mpmath as mp
import numpy as np
from aoki_majorana import uv, greater, contraction, slash, I, EPSILON_LOWER as EPS, J
from weyl_scalar_seed import scalar_seed_folded
mp.mp.dps=40
lines=[]

def emit(text):
    lines.append(text);print(text,flush=True)


def old_modes(mu,x,h):
    low=mu/mp.sqrt(2*x)*mp.whitw(-mp.mpf('.5'),1j*mu,-2j*x)
    high=mp.whitw(mp.mpf('.5'),1j*mu,-2j*x)/mp.sqrt(2*x)
    return (low,high) if h==1 else (high,low)


def eom_residual(modes,mu,x,h):
    u,v=modes(mu,x,h)
    first=-1j*mp.diff(lambda y:modes(mu,y,h)[0],x)+h*u-mu/x*v
    second=-1j*mp.diff(lambda y:modes(mu,y,h)[1],x)-h*v-mu/x*u
    return max(abs(first),abs(second))

for mu in map(mp.mpf,('.2','.73','2')):
    current=[];old=[]
    for x,h in product(map(mp.mpf,('.3','1','10')),(1,-1)):
        current.append(eom_residual(uv,mu,x,h))
        old.append(eom_residual(old_modes,mu,x,h))
    assert max(current)<mp.mpf('1e-35')
    assert max(old)>mp.mpf('.1')
    uplus,high=old_modes(mu,mp.mpf('100'),1)
    ratio=abs(uplus/high)
    assert ratio<mp.mpf('.02')
    emit(f'CHECK modes mu={mu}: max EOM residual v1={float(max(old)):.6g}, v2/Aoki={float(max(current)):.3e}; |u+/u-| at x=100 is {float(ratio):.8f}, not 1 as v7_24 (13) requires')

n=np.array([.3,-.4,np.sqrt(.75)])
for mu in (.2,.73,2.):
    x=.7;k=n;t=-x
    gt=greater(k,t,t,mu)[:2,2:]
    wrong=np.zeros((2,2),complex)
    for h in (1,-1):
        u,v=map(complex,uv(mu,x,h))
        wrong-=abs(v)**2*(I-h*slash(n))/2
    correct=-greater(-k,t,t,mu)[2:,:2].T
    bad=np.linalg.norm(gt-wrong-I)
    good=np.linalg.norm(gt-correct-I)
    assert good<1e-13 and bad>.05
    emit(f'CHECK canonical normal anticommutator mu={mu}, x=.7: missing -k residual={bad:.8f}, corrected residual={good:.3e}')

for mu in map(mp.mpf,('.2','.73','2','5')):
    # You v2 (2.39a) at zero chemical potential, compare the coefficient
    # of (k^2 tau1 tau2)^(+i mu) to Aoki's c/2 via Gamma duplication.
    b=mp.gamma(-2j*mu)**2/mp.gamma(-1j*mu)**2*4**(1j*mu)
    c=mp.power(2,-1-2j*mu)*mp.gamma(.5-1j*mu)**2/mp.pi
    rel=abs(b-c/2)/abs(b)
    assert rel<mp.mpf('1e-35')
    gt=-complex(b)*slash(n)
    # P_{-s}(-k)=P_s(k) preserves the soft normal contraction.
    lt=gt.copy()
    wrong_lt=-gt
    good_s=(gt+lt)/2;bad_s=(gt+wrong_lt)/2
    assert np.linalg.norm(bad_s)==0 and np.linalg.norm(good_s)>0
    emit(f'CHECK soft normal mu={mu}: You-v2/Aoki coefficient relative difference={float(rel):.3e}; ||D_S,NL,+||/(k^2 tau1 tau2)^i mu magnitude={np.linalg.norm(good_s):.8g}; v7_24 (18) gives 0')

left0=np.vstack((I,EPS));right0=np.hstack((EPS,I))
q=np.array([.2,-.3,.4]);K=q+np.array([0.,0.,1.])
sq=slash(q/np.linalg.norm(q));sK=slash(K/np.linalg.norm(K))
nL=np.array([.6,0.,.8]);nR=np.array([.3,np.sqrt(.75),-.4])
for mu in (.2,.73,2.):
    residual=0.;single_max=0.
    for a,b,c,d in product((1,-1),repeat=4):
        vals=[]
        for signL,signR in product((1,-1),repeat=2):
            DL=contraction(signL*nL,-.37,-1.41,a,b,mu)
            DR=contraction(signR*nR,-.61,-1.17,c,d,mu)
            VL=-.5*right0@J@DL@J@left0
            VR=-.5*right0@J@DR@J@left0
            vals.append(np.trace(VL@sK@VR@sq))
        residual=max(residual,abs(sum(vals)))
        single_max=max(single_max,abs(vals[0]))
    assert residual<1e-13 and single_max>.001
    emit(f'CHECK Bose cancellation mu={mu}, all 16 SK: largest single-cycle spin coefficient={single_max:.8g}, largest four-ordering sum={residual:.3e}')

# A reproducible branch selection diagnostic for the SAME physical hard kernel.
# Every branch contains the same endpoint regulator delta. No branch is assigned
# its own finite part. The normalization is H_GF=sum_ab H_GF^(ab), including ab.
mp.mp.dps=25
mu=mp.mpf('.73');nu=mu-.5j
for delta in map(mp.mpf,('.8','.2','.1')):
    alpha=delta+1j*mu
    terms={}
    for a,b in product((1,-1),repeat=2):
        terms[a,b]=-sum((1j*a)**j*(1j*b)**l*
                       scalar_seed_folded(a,b,alpha-2+j,alpha-2+l,nu)
                       for j,l in product((0,1),repeat=2))
    full=sum(terms.values())
    emit(f'SK example mu=.73 delta={delta}: ++={mp.nstr(terms[1,1],11)}, +-={mp.nstr(terms[1,-1],11)}, -+={mp.nstr(terms[-1,1],11)}, --={mp.nstr(terms[-1,-1],11)}; full={mp.nstr(full,11)}')
    assert abs(terms[1,1]-full)>.1

emit('Scope: printed mode/EOM and canonical counterexamples; nonzero soft normal coefficient; original-field zeroth-order exchange cancellation. SK table reuses the scalar evaluator. No complete audit of You v2 or numerical reconstruction of v7_24 4F3 coefficients.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
