"""Independent algebraic/numerical checks of Weyl Spinor Signal.tex.
Run: python3 checks/weyl_signal_checks.py
No claim of numerical continuation of the full hard SK seeds is made.
"""
import mpmath as mp
import numpy as np
import sympy as sp
from pathlib import Path
mp.mp.dps = 40
report=[]
def check(name, err, tol=1e-10):
    err=float(err)
    assert err < tol, (name,err,tol)
    report.append(f'PASS {name}: residual {err:.3e}')
mu=mp.mpf('.73')
F=lambda x:mp.sqrt(mp.pi*x)/2*mp.exp(mp.pi*mu/2)*mp.hankel1(.5-1j*mu,x)
G=lambda x:mp.sqrt(mp.pi*x)/2*mp.exp(-mp.pi*mu/2)*mp.hankel1(.5+1j*mu,x)
errs=[]
for x in map(mp.mpf,['.001','.7','12']):
    errs += [abs(-1j*mp.diff(F,x)-mu/x*F(x)-G(x)),abs(-1j*mp.diff(G,x)+mu/x*G(x)-F(x)),abs(abs(F(x))**2+abs(G(x))**2-1)]
check('mode equations and normalization',max(errs),1e-30)
A=-1j*mp.exp(mp.pi*mu/2)*2**(-.5-1j*mu)*mp.gamma(.5-1j*mu)/mp.sqrt(mp.pi)
B=-1j*mp.exp(-mp.pi*mu/2)*2**(-.5+1j*mu)*mp.gamma(.5+1j*mu)/mp.sqrt(mp.pi)
d=1j/(1+2j*mu); dm=1j/(1-2j*mu)
for x in [mp.mpf('1e-7')]:
    check('F first descendant',abs(F(x)-A*x**(1j*mu)-dm*B*x**(1-1j*mu)),1e-13)
    check('G first descendant',abs(G(x)-B*x**(-1j*mu)-d*A*x**(1+1j*mu)),1e-13)

def mode_parts(x,nu,N):
    return N*mp.sqrt(x)*mp.besselj(-nu,x)/(1j*mp.sin(mp.pi*nu)), -N*mp.sqrt(x)*mp.exp(-1j*mp.pi*nu)*mp.besselj(nu,x)/(1j*mp.sin(mp.pi*nu))
def parts(x):
    fp,fm=mode_parts(x,.5-1j*mu,mp.sqrt(mp.pi)/2*mp.exp(mp.pi*mu/2))
    gm,gp=mode_parts(x,.5+1j*mu,mp.sqrt(mp.pi)/2*mp.exp(-mp.pi*mu/2))
    return fp,fm,gp,gm
x,y=mp.mpf('.27'),mp.mpf('.49')
fp,fm,gp,gm=parts(x); f2p,f2m,g2p,g2m=parts(y)
up=mp.matrix([fp,gp]); um=mp.matrix([fm,gm]); u2p=mp.matrix([f2p,g2p]);u2m=mp.matrix([f2m,g2m])
vp=mp.matrix([mp.conj(gm),-mp.conj(fm)]); vm=mp.matrix([mp.conj(gp),-mp.conj(fp)])
v2p=mp.matrix([mp.conj(g2m),-mp.conj(f2m)]);v2m=mp.matrix([mp.conj(g2p),-mp.conj(f2p)])
check('exact nonanalytic Wightman equality including descendants',mp.norm(up*u2m.H+um*u2p.H+vp*v2m.H+vm*v2p.H),1e-30)
u=mp.matrix([F(x),G(x)]);v=mp.matrix([mp.conj(G(x)),-mp.conj(F(x))])
check('signed equal-time completeness',mp.norm(u*u.H+v*v.H-mp.eye(2)),1e-30)

# Direct Pauli contraction and transverse angular average, independent of Gamma continuation.
sigma=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
def slash(v):return np.einsum('i,ijk->jk',v,sigma)
rng=np.random.default_rng(730)
worst=0
for _ in range(20):
    aa,bb,q,K=rng.normal(size=(4,3));q/=np.linalg.norm(q);K/=np.linalg.norm(K)
    lhs=np.trace(slash(aa)@slash(K)@slash(bb)@slash(q))
    rhs=2*((aa@K)*(bb@q)+(aa@q)*(bb@K)-(aa@bb)*(q@K))
    worst=max(worst,abs(lhs-rhs))
check('four-Pauli trace',worst)
a=.7+.2j;b=-.3+.8j
n=np.array([np.sin(.7),0,np.cos(.7)]);m=np.array([np.sin(1.1)*np.cos(.9),np.sin(1.1)*np.sin(.9),np.cos(1.1)])
c1=n[2];c3=m[2];e=n@m
rr=.61;w=.37;L=rr+w*w
M=L*L+L/4-w*w/2;P=(L+.25)*rr/2;Q=(L+.25)*(w*w-rr/2)-w*w/2
E=rr/2;FF=w*w-rr/2;T=rr*rr/8;U=rr*w*w/2-T;V=w**4-6*U-3*T
X=(L-.25)*rr/2;Y=(L-.25)*(w*w-rr/2)
rhs=a*a*M+a*b*(2*P+Q*(c1*c1+c3*c3))+b*b*(2*T+(4*T-X)*e*e+2*U*(c1*c1+c3*c3)+(8*U-Y-E/2)*e*c1*c3+(2*V-FF/2)*c1*c1*c3*c3)
sumv=0j
for phi in np.arange(512)*2*np.pi/512:
    ell=np.array([np.sqrt(rr)*np.cos(phi),np.sqrt(rr)*np.sin(phi),w]);z=np.array([0,0,1])
    aa=a*ell+b*n*(n@ell);bb=a*ell+b*m*(m@ell)
    sumv+=2*(aa@ell)*(bb@ell)-.5*(aa@z)*(bb@z)-(aa@bb)*(L-.25)
check('azimuthal bubble polynomial',abs(sumv/512-rhs))
w0=a*a*M+2*a*b*P+2*b*b*T;w2=a*b*Q+2*b*b*U;we=b*b*(4*T-X);wez=b*b*(8*U-Y-E/2);w4=b*b*(2*V-FF/2)
ang=(w0+we/2)+(w2-we/2)*(c1*c1+c3*c3)+(w4+wez+1.5*we)*c1*c1*c3*c3+(2*we+wez)/4*np.sin(1.4)*np.sin(2.2)*np.cos(.9)+we/2*np.sin(.7)**2*np.sin(1.1)**2*np.cos(1.8)
check('five angular coefficients',abs(ang-rhs))
# Bose-odd zeroth-order block must vanish at either side.
q=np.array([.2,.3,-.1]);K=q+np.array([0,0,1]);q/=np.linalg.norm(q);K/=np.linalg.norm(K)
check('zeroth-order four-cycle cancellation',abs(sum(np.trace(slash(sn*n)@slash(K)@slash(sm*m)@slash(q)) for sn in [-1,1] for sm in [-1,1])))

# Independent check of ordered Laplace kernel after x=z*y.
alpha=mp.mpf('1.27');beta=mp.mpf('.91');s=mp.mpf('2.3');t=mp.mpf('1.7')
val=mp.gamma(alpha+beta)*mp.quad(lambda z:z**(alpha-1)/(t+s*z)**(alpha+beta),[0,1])
closed=mp.gamma(alpha+beta)/(alpha*t**(alpha+beta))*mp.hyp2f1(alpha+beta,alpha,alpha+1,-s/t)
check('ordered scalar Laplace integral',abs(val-closed),1e-30)
# Validate coefficient normalization of both Hankel-mode series.
for nu,N,func in [(.5-1j*mu,mp.sqrt(mp.pi)/2*mp.exp(mp.pi*mu/2),F),(.5+1j*mu,mp.sqrt(mp.pi)/2*mp.exp(-mp.pi*mu/2),G)]:
    xx=mp.mpf('.6');val=0
    for eta in [-1,1]:
        for j in range(30):
            pref=N/(1j*mp.sin(mp.pi*nu)) if eta==-1 else -N*mp.exp(-1j*mp.pi*nu)/(1j*mp.sin(mp.pi*nu))
            val+=pref*(-1)**j*2**(-eta*nu-2*j)/(mp.factorial(j)*mp.gamma(j+1+eta*nu))*xx**(.5+eta*nu+2*j)
    check('Hankel power-series normalization',abs(val-func(xx)),1e-30)
# Output finite monomials for the bubble functional J.
u,v=sp.symbols('u v');ww=(v-u)/2;LL=(u+v)/2-sp.Rational(1,4);RR=LL-ww**2
polys={'M':LL**2+LL/4-ww**2/2,'P':(LL+sp.Rational(1,4))*RR/2,'Q':(LL+sp.Rational(1,4))*(ww**2-RR/2)-ww**2/2,'E':RR/2,'F':ww**2-RR/2,'T':RR**2/8}
polys['U']=RR*ww**2/2-polys['T'];polys['V']=ww**4-6*polys['U']-3*polys['T'];polys['X']=(LL-sp.Rational(1,4))*RR/2;polys['Y']=(LL-sp.Rational(1,4))*(ww**2-RR/2)
# Compare independently evaluated Gamma sums with the compact rational moments.
xx=mp.mpf('.5')-1j*mu
f0=lambda aa,bb:mp.gamma(aa+bb-mp.mpf('1.5'))*mp.gamma(mp.mpf('1.5')-aa)*mp.gamma(mp.mpf('1.5')-bb)/(mp.gamma(aa)*mp.gamma(bb)*mp.gamma(3-aa-bb))
base=f0(xx,xx)
rat={
'M':3*(xx-1)*(2*xx-3)/(8*(xx-2)*(4*xx-7)*(4*xx-5)),
'P':(xx-1)*(2*xx-3)/(8*(xx-2)*(4*xx-7)*(4*xx-5)),
'Q':0,
'E':(2*xx-3)/(8*(xx-2)*(4*xx-5)),
'F':-1/(8*(xx-2)),
'T':(2*xx-5)*(2*xx-3)/(64*(xx-3)*(xx-2)*(4*xx-7)*(4*xx-5)),
'U':-(2*xx-3)/(64*(xx-3)*(xx-2)*(4*xx-5)),
'V':3/(64*(xx-3)*(xx-2)),
'X':-(2*xx-5)*(2*xx-3)/(16*(xx-2)*(4*xx-7)*(4*xx-5)),
'Y':(2*xx-3)/(16*(xx-2)*(4*xx-5))}
for name,poly in polys.items():
    val=sum(mp.mpf(str(coeff))*f0(xx-m,xx-n) for (m,n),coeff in sp.Poly(poly,u,v).terms())
    check('Gamma recurrence for '+name,abs(val-base*rat[name]),1e-30)
aa=mp.mpf('.5')+1j*mu;bb=mp.mpf('.5')-1j*mu
jm=(f0(aa-1,bb)+f0(aa,bb-1)-f0(aa,bb))/2
check('mixed-branch Gamma bubble',abs(jm-2*mp.sqrt(mp.pi)/3*mu*(1+mu**2)*mp.coth(mp.pi*mu)),1e-30)
report.append('\nBubble moments: replace u**m*v**n with F0(1/2-i*mu-m,1/2-i*mu-n).')
for name,poly in polys.items(): report.append(f'{name} = {sp.expand(poly)}')
report.append('\nScope: algebra and convergent scalar identities checked; full hard-seed analytic continuation has not been numerically scanned.')
out='\n'.join(report)+'\n'
Path(__file__).with_name('weyl_signal_checks.txt').write_text(out)
print(out)
