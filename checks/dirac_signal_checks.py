"""Verify the four-component soft expansion and scalar integrals.

Run: PYTHONDONTWRITEBYTECODE=1 python3 checks/dirac_signal_checks.py
Independent mode sums test the scalar-kernel projections and the four-cut
soft expansion. Scalar convolution checks reuse the existing Gamma engine;
new source parsing tests the actual J expressions. Time checks compare the
single-3F2 finite form with the Q family and with complete regulated SK sums.
This is a leading collapsed same-branch check, not a full box integration.
"""
from contextlib import redirect_stdout
from functools import lru_cache
from io import StringIO
from itertools import product
from pathlib import Path
import re
import mpmath as mp
import numpy as np
import sympy as s
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                       implicit_multiplication_application)
from dirac_propagator_checks import (IDENTITY, G0, GAMMA, spatial_gamma,
                                     scalar_kernel, select_sk)
from dirac_box_cut_checks import full_line, cut_line, endpoint, soft_prefactor
import top_loop_j_checks as loop
from whittaker_scalar_seed import hard_coefficients, extrapolated_hard_coefficients

ROOT = Path(__file__).resolve().parents[1]


def scalar_f(p,t1,t2,mu,eta,a,b):
    """Apply the displayed scalar differential operators, independently of modes."""
    p=mp.mpf(p); X,Y=-p*t1,-p*t2
    def kernel(x,kap):
        return select_sk(scalar_kernel(p,x,Y,mu,kap,False),
                         scalar_kernel(p,x,Y,mu,kap,True),t1,t2,a,b)
    dm,dp=kernel(X,-.5),kernel(X,.5)
    rm=-X*mp.diff(lambda x:kernel(x,-.5),X)-1j*X*dm+dm/2
    ap=X*mp.diff(lambda x:kernel(x,.5),X)-1j*X*dp-dp/2
    z=1j*eta*mu
    return tuple(map(complex,(-1j*p/mp.sqrt(X*Y)*(ap-z*dp+rm+z*dm),
                             -1j*p/mp.sqrt(X*Y)*(rm-z*dm-ap-z*dp))))


def mode_f(p,n,t1,t2,mu,eta,a,b):
    d=full_line(tuple(p*n),t1,t2,mu,a,b)
    lp,lm=(IDENTITY+eta*G0)/2,(IDENTITY-eta*G0)/2
    return -np.trace(lm@d@spatial_gamma(n)),-eta*np.trace(lp@d)


def hard_v(e,k,n,q,K,t1,t2,mu,eta,a,b):
    ell=(q+K)/2
    en=lambda p,t: IDENTITY if e==0 else endpoint(e*p,t,mu,eta)
    return (en(q,t1)@full_line(tuple(k*n+e*ell),t1,t2,mu,a,b)@en(K,t2)
           +en(q,t2)@full_line(tuple(-k*n+e*ell),t2,t1,mu,b,a)@en(K,t1))


def velocity(k,n,q,K,t1,t2,mu,eta,a,b):
    h=2e-5*k
    f0,f1=mode_f(k,n,t1,t2,mu,eta,a,b)
    fr=mode_f(k,n,t2,t1,mu,eta,b,a)[1]
    df=(mode_f(k+h,n,t1,t2,mu,eta,a,b)[0]
       -mode_f(k-h,n,t1,t2,mu,eta,a,b)[0])/(2*h)
    pi1=f0-1j*k/(1+2j*eta*mu)*(t1*f1+t2*fr)
    pi2=k*df-f0
    ell=(q+K)/2
    return pi1*ell+pi2*n*np.dot(n,ell)


def check_projection_and_trace():
    gm=[s.Matrix(g).applyfunc(s.nsimplify) for g in GAMMA]
    for eta in (1,-1):
        lm=(s.eye(4)-eta*gm[0])/2
        for i,j,k,l in product(range(1,4),repeat=4):
            want=2*(int(i==j)*int(k==l)-int(i==k)*int(j==l)+int(i==l)*int(j==k))
            assert s.trace(lm*gm[i]*gm[j]*gm[k]*gm[l])==want
    print('PASS exact projected four-gamma identity for both soft branches')
    n=np.array([.6,0,.8]); k=1.3; q=np.array([.3,-.2,.1]);K=q+[0,0,.7]
    maxima=np.zeros(4); raw=0
    for mu,eta,a,b in product((.2,.73,2.),(1,-1),(1,-1),(1,-1)):
        for t1,t2 in ((-.43,-1.27),(-1.27,-.43)):
            lp,lm=(IDENTITY+eta*G0)/2,(IDENTITY-eta*G0)/2
            f=mode_f(k,n,t1,t2,mu,eta,a,b)
            fs=scalar_f(k,t1,t2,mu,eta,a,b)
            maxima[0]=max(maxima[0],np.linalg.norm(np.array(f)-fs))
            rev=mode_f(k,n,t2,t1,mu,eta,b,a)[0]
            maxima[1]=max(maxima[1],abs(f[0]-rev))
            v0=hard_v(0,k,n,q,K,t1,t2,mu,eta,a,b)
            raw=max(raw,np.linalg.norm(v0))
            maxima[2]=max(maxima[2],np.linalg.norm(lm@v0@lp))
            h=2e-5
            dv=(hard_v(h,k,n,q,K,t1,t2,mu,eta,a,b)
               -hard_v(-h,k,n,q,K,t1,t2,mu,eta,a,b))/(2*h)
            v=velocity(k,n,q,K,t1,t2,mu,eta,a,b)
            maxima[3]=max(maxima[3],np.linalg.norm(lm@dv@lp-lm@spatial_gamma(v)@lp/k))
    assert max(maxima[:3])<2e-13 and maxima[3]<2e-8, maxima
    print('PASS modes vs scalar F0,F1; reversal; projected zero; first derivative')
    print('  3 masses, both branches, all four hard SK choices, both time orders')
    print('  maximum residuals:',*[f'{v:.3g}' for v in maxima])
    print(f'  unprojected zeroth-order sum is nonzero (sample maximum {raw:.6g})')


def check_full_cut_limit():
    q=np.array([.3,-.2,.1]);K=q+[0,0,.7];ell=(q+K)/2;ks=K-q
    nL=np.array([.6,0,.8]); nR=np.array([.3,.4,np.sqrt(.75)])
    kl,kr=1.3,1.7; ts=[-.43,-1.27,-.71,-1.59]
    maxerr=0.; ratios=[]
    for mu,eta in product((.2,.73),(1,-1)):
        vl={};vr={}
        for a,b in product((1,-1),repeat=2):
            vl[a,b]=velocity(kl,nL,q,K,*ts[:2],mu,eta,a,b)
            vr[a,b]=velocity(kr,nR,K,q,*ts[2:],mu,eta,a,b)
        for aa in product((1,-1),repeat=4):
            vL,vR=vl[aa[:2]],vr[aa[2:]]
            numerator=(2*np.dot(vL,ell)*np.dot(vR,ell)-.5*np.dot(vL,ks)*np.dot(vR,ks)
                      -np.dot(vL,vR)*(np.dot(ell,ell)-np.dot(ks,ks)/4))
            leading=8*numerator/(kl*kr*np.linalg.norm(q)*np.linalg.norm(K))
            errors=[]
            # Each projected block starts at O(e), with a generic O(e^2)
            # remainder. Their product therefore has O(e^3) corrections:
            # the relative error of the leading scalar trace is O(e).
            for e in (.01,.005,.0025):
                def D(i,j,p):return full_line(tuple(p),ts[i-1],ts[j-1],mu,aa[i-1],aa[j-1])
                def S(i,j,p):return cut_line(tuple(p),ts[i-1],ts[j-1],mu,eta)
                left,right=kl*nL+e*ell,kr*nR+e*ell
                lrev,rrev=-kl*nL+e*ell,-kr*nR+e*ell
                qq,kk=e*q,e*K
                actual=np.trace(D(1,2,left)@S(2,3,kk)@D(3,4,right)@S(4,1,qq)
                    +D(1,2,left)@S(2,4,kk)@D(4,3,rrev)@S(3,1,qq)
                    +D(2,1,lrev)@S(1,3,kk)@D(3,4,right)@S(4,2,qq)
                    +D(2,1,lrev)@S(1,4,kk)@D(4,3,rrev)@S(3,2,qq))
                pref=soft_prefactor(kk,ts[1],ts[2],mu,eta)*soft_prefactor(qq,ts[3],ts[0],mu,eta)
                err=abs(actual/pref/e**2-leading)/max(abs(leading),1e-8)
                errors.append(err)
            maxerr=max(maxerr,errors[-1]);ratios.extend([errors[0]/errors[1],errors[1]/errors[2]])
    assert maxerr<.003 and min(ratios)>1.8 and max(ratios)<2.2,(maxerr,min(ratios),max(ratios))
    print('PASS exact four physical cut traces -> displayed scalar numerator')
    print('  2 masses, both branches, all 16 SK assignments; soft scales .01,.005,.0025')
    print(f'  max relative residual at .0025: {maxerr:.3g}; error ratios {min(ratios):.5f}--{max(ratios):.5f}')


def check_branch_conjugation():
    n=np.array([.6,0,.8]);k=1.3;t1,t2=-.43,-1.27
    q=np.array([.3,-.2,.1]);K=q+[0,0,.7]
    worst=0.
    for mu,a,b in product((.2,.73,2.),(1,-1),(1,-1)):
        fm=np.array(mode_f(k,n,t1,t2,mu,-1,a,b))
        fp=np.conj(mode_f(k,n,t1,t2,mu,1,-a,-b))
        worst=max(worst,np.linalg.norm(fm-fp*np.array([1,-1])))
        vm=velocity(k,n,q,K,t1,t2,mu,-1,a,b)
        vp=velocity(k,n,q,K,t1,t2,mu,1,-a,-b)
        worst=max(worst,np.linalg.norm(vm-vp.conj()))
    assert worst<1e-9,worst
    print(f'PASS soft-branch conjugation at reversed SK signs, residual {worst:.3g}')


def check_single_ordering():
    n=np.array([.6,0,.8]);k=1.3;t1,t2=-.43,-1.27
    q=np.array([.3,-.2,.1]);K=q+[0,0,.7];ell=(q+K)/2
    worst=0.
    for mu,eta,a,b in product((.2,.73,2.),(1,-1),(1,-1),(1,-1)):
        lm,lp=(IDENTITY-eta*G0)/2,(IDENTITY+eta*G0)/2
        f0,f1=mode_f(k,n,t1,t2,mu,eta,a,b)
        fr=mode_f(k,n,t2,t1,mu,eta,b,a)[1]
        h=2e-5
        df=(mode_f(k+h,n,t1,t2,mu,eta,a,b)[0]
           -mode_f(k-h,n,t1,t2,mu,eta,a,b)[0])/(2*h)
        predicted=(f0*ell+(k*df-f0)*n*np.dot(n,ell)
                  -1j*k/(1+2j*eta*mu)*(t1*f1*q+t2*fr*K))/(2*k)
        def single(e):
            return (endpoint(e*q,t1,mu,eta)
                    @full_line(tuple(k*n+e*ell),t1,t2,mu,a,b)
                    @endpoint(e*K,t2,mu,eta))
        actual=lm@(single(h)-single(-h))@lp/(2*h)
        worst=max(worst,np.linalg.norm(actual-lm@spatial_gamma(predicted)@lp))
    assert worst<1e-8,worst
    print(f'PASS single-ordering hard and endpoint correction formulas, residual {worst:.3g}')


def check_source_j():
    tex=(ROOT/'Top Quark Signal.tex').read_text()
    end=tex.index(r'\label{eq:dirac-loop-reduction}')
    block=tex[tex.rfind(r'\begin{align}',0,end):end]
    lhs12=r'\mathcal I^{12}_{t}';lhs22=r'\mathcal I^{22}_{t}'
    start=block.index(lhs12,block.index(r'\nonumber\\'))
    mid=block.index(lhs22)
    terms={'P':block[start+len(lhs12):mid],'H':block[mid+len(lhs22):]}
    variables={}
    angular={r'(\cos^2\theta_1+\cos^2\theta_3-3\cos^2\theta_1\cos^2\theta_3)':'angleA',
             r'\sin(2\theta_1)\sin(2\theta_3)\cos\phi_3':'angleB',
             r'\sin^2\theta_1\sin^2\theta_3\cos(2\phi_3)':'angleD'}
    variables.update({v:s.Symbol(v) for v in angular.values()})
    def symbol(m):
        r,p,n=map(int,m.groups());assert n==1
        name=f'j{r}m{p}n{n}';variables[name]=s.Symbol(name);return '('+name+')'
    def formal(poly):
        return sum(c*s.Symbol(f'j{r}m{p-1}n1') for (p,r),c in s.Poly(poly,loop.q,loop.z).terms())
    expect={'P':formal(loop.polys['P']), 'H':formal(loop.polys['H0'])
            +variables['angleA']*formal(loop.polys['H2'])
            +variables['angleB']*formal(loop.polys['HB'])
            +variables['angleD']*formal(loop.polys['HD'])}
    for name,rhs in terms.items():
        for old in ('&=',r'\nonumber\\','&',r'\quad',r'\Big'):rhs=rhs.replace(old,'')
        rhs=re.sub(r'\\hspace\{[^}]+\}','',rhs)
        rhs=' '.join(rhs.split())
        for old,new in angular.items():rhs=rhs.replace(old,'('+new+')')
        rhs=re.sub(r'\\mathcal J_(\d)\((\d+),(\d+)\)',symbol,rhs)
        rhs=re.sub(r'\\frac(\{[0-9]+\}|[0-9])(\{[0-9]+\}|[0-9])',
                   lambda m:f"({m[1].strip('{}')}/{m[2].strip('{}')})",rhs)
        rhs=rhs.replace('[','(').replace(']',')').strip().rstrip(',.')
        result=parse_expr(rhs,local_dict=variables,transformations=
                          standard_transformations+(implicit_multiplication_application,))
        assert s.expand(result-expect[name])==0
    print('PASS actual Appendix A J sums against azimuthal numerator polynomials')


def finite_p(nu):
    z=1j*nu
    f=mp.hyper([1,1+z,1+2*z],[2+2*z,2+2*z],1)
    e=mp.pi**2*mp.gamma(1+2*z)/(mp.cosh(mp.pi*nu)*mp.gamma(1+z)**2)
    g=mp.gamma(1+2*z)*mp.cosh(mp.pi*nu)
    return (g*(e-2*f/(1+2*z)**2-(5+2*z)*(4+9*z+4*z*z)/(2*(1+z)*(2+z)*(1+2*z))),
            g*(2*f/(1+2*z)-(1+2*z)*e+z*(5+2*z)/(2*(1+z)*(2+z))))


def check_scalar_source_reuse():
    new=(ROOT/'Top Quark Signal.tex').read_text()
    sms=(ROOT/'SMS in CHC.tex').read_text()
    def equation(tex,prefix,name):
        tex=re.sub(r'(?m)^\s*%.*$', '',tex)
        end=tex.index(r'\label{eq:'+prefix+name+'}')
        start=max(tex.rfind(r'\begin{align}',0,end),
                  tex.rfind(r'\begin{equation}',0,end))
        out=tex[start:end].replace('eq:top-','eq:dirac-')
        out=out.replace(r'\tn',r'\tilde\nu')
        out=re.sub(r'\\Dsc\{([^}]+)\}',r'D^{(\1)}',out)
        return re.sub(r'\s+','',out)
    labels=('hard-kernels','loop-reduction','loop-Gamma','loop-closed',
            'K-definition','K-scaling','K-opposite','K-same','H-definition',
            'H0-K','H1-K','K-derivative','P-from-seed','finite-master',
            'P-explicit-finite','common-factor','angular-closed')
    for name in labels:
        assert equation(new,'dirac-',name)==equation(sms,'top-',name),name
    print(f'PASS {len(labels)} scalar equation blocks identical to SMS after label/macro normalization')


def check_time_and_angles():
    mp.mp.dps=30
    worst=mp.mpf(0)
    for nu in map(mp.mpf,('.2','.73','2')):
        new=finite_p(nu);old=hard_coefficients(nu,30)
        # The pre-existing scalar helper has F_j(old)=-F_j(SMS)/2.
        # This follows directly by comparing its differential operators;
        # 2^(2 i nu) then converts k_L to k_12, as in the new appendix.
        mapped=[-2*mp.power(2,2j*nu)*v for v in old[:2]]
        err=max(abs(a-b)/max(1,abs(a)) for a,b in zip(new,mapped))
        worst=max(worst,err);assert err<mp.mpf('1e-25')
    print('PASS finite single-3F2 coefficients vs independent Q-family reduction')
    print('  3 masses; relative residual',mp.nstr(worst,5))
    nu=mp.mpf('.73');new=finite_p(nu)
    regulated=extrapolated_hard_coefficients(nu,30,step='.0008',levels=4)
    mapped=[-2*mp.power(2,2j*nu)*v for v in regulated[:2]]
    err=max(abs(a-b)/max(1,abs(a)) for a,b in zip(new,mapped))
    assert err<mp.mpf('1e-9'),err
    print('PASS complete regulated SK/kappa sums vs final finite P1,P2 at nu=.73')
    print('  four common regulators; relative residual',mp.nstr(err,5))
    z,p1,p2=s.symbols('z p1 p2')
    # Divide all angular coefficients by the common B_t. Gamma factors cancel.
    common=(5+2*z)*(3+2*z)*(5+4*z)*(3+4*z)
    i12=-(1+z)*(1+2*z)/(4*(3+2*z)*(5+4*z)*(3+4*z))
    h=-(1+z)/(4*common)
    vector=[8*common*((3*p1**2+2*p1*p2)*i12+p2**2*h*(2+z)*(1+2*z)),
            8*common*p2**2*h*(-(1+z)*(1+2*z)),
            8*common*p2**2*h*(3*(1+z)*(1+2*z)),
            8*common*p2**2*h*(-(1+2*z)*(3+2*z)/2),
            8*common*p2**2*h*((2+z)*(3+2*z))]
    want=[-2*(1+z)*(1+2*z)*((5+2*z)*(3*p1**2+2*p1*p2)+(2+z)*p2**2),
          2*p2**2*(1+z)**2*(1+2*z),-6*p2**2*(1+z)**2*(1+2*z),
          p2**2*(1+z)*(1+2*z)*(3+2*z),-2*p2**2*(1+z)*(2+z)*(3+2*z)]
    assert all(s.factor(a-b)==0 for a,b in zip(vector,want))
    # Four cut traces (8), two time integrals (1/16), 1/(kL kR) -> 4/(k12 k34).
    assert -s.Rational(8,16)*4 == -s.Rational(1,4)*8
    print('PASS five angular coefficients, A3=-3 A2, final overall normalization')


def main():
    mp.mp.dps=30
    log=StringIO()
    class Tee:
        def write(self,text):
            import sys
            sys.__stdout__.write(text);sys.__stdout__.flush();log.write(text)
        def flush(self):pass
    with redirect_stdout(Tee()):
        check_projection_and_trace()
        check_full_cut_limit()
        check_branch_conjugation()
        check_single_ordering()
        loop.check_azimuth();check_source_j()
        loop.check_reduction();loop.check_nu_form()
        check_scalar_source_reuse()
        check_time_and_angles()
        print('Scope: leading same-branch collapsed signal; no direct unexpanded box integral.')
    Path(__file__).with_suffix('.txt').write_text(log.getvalue())

if __name__=='__main__':main()
