"""Independent BD hard-time quadrature at positive common endpoint regulator.
No hypergeometric scalar seed enters quadrature. Same-SK ordered triangles
are Wick rotated together; opposite-SK integrals use separate BD contours.
This tests physical fermion mass indices, but not the delta->0 limit.
"""
from pathlib import Path
import numpy as np
import mpmath as mp
from weyl_scalar_seed import scalar_seed_folded, regulated_diagonal_seed, regulated_hard_coefficients
mp.mp.dps=28


def modes(mu,x):
    n=mp.sqrt(mp.pi*x)/2
    return np.array([complex(n*mp.exp(mp.pi*mu/2)*mp.hankel1(.5-1j*mu,x)),
                     complex(n*mp.exp(-mp.pi*mu/2)*mp.hankel1(.5+1j*mu,x)),
                     complex(n*mp.exp(mp.pi*mu/2)*mp.hankel2(.5+1j*mu,x)),
                     complex(n*mp.exp(-mp.pi*mu/2)*mp.hankel2(.5-1j*mu,x))])


def integrate(mu,alpha,beta,n=40,limit=22.,diagonal_shift=0):
    nodes,w=np.polynomial.legendre.leggauss(n)
    # Cubic endpoint maps smooth the integrable regulated powers.
    u=(nodes+1)/2*limit**(1/3); wu=w/2*limit**(1/3)
    ts=u**3;wt=wu*3*u**2
    v=(nodes+1)/2;wv=w/2
    answer=np.zeros(3,complex)
    lap={}; dlap={}
    def radial_derivative(values,arguments):
        f,g,fb,gb=values.T
        return np.array([1j*mu*f+1j*arguments*g,
                         1j*arguments*f-1j*mu*g,
                         -1j*mu*fb-1j*arguments*gb,
                         -1j*arguments*fb+1j*mu*gb]).T
    for sign in (1,-1):
        rotation=-1j*sign
        outer=np.array([modes(mu,rotation*t) for t in ts])
        # E_sign(rotation*t)=(1+t)*exp(-t).
        external=(1+ts)*np.exp(-ts)
        douter=radial_derivative(outer,rotation*ts)
        for power in (alpha,beta,alpha+diagonal_shift):
            lap[sign,power]=rotation**power*np.sum((wt*ts**(power-1)*external)[:,None]*outer,axis=0)
            dlap[sign,power]=rotation**power*np.sum((wt*ts**(power-1)*external)[:,None]*douter,axis=0)
        for t,weight,large,dlarge in zip(ts,wt,outer,douter):
            ss=t*v**3;ws=t*3*v**2*wv
            small=np.array([modes(mu,rotation*s) for s in ss])
            dsmall=radial_derivative(small,rotation*ss)
            w1=ss**(alpha-1)*t**(beta-1)
            w2=ss**(beta-1)*t**(alpha-1)
            ext=(1+ss)*np.exp(-ss)*(1+t)*np.exp(-t)
            if sign==1:
                h21=(w1+w2)*small[:,1]*large[2]
                dh21=(w1+w2)*(dsmall[:,1]*large[2]+small[:,1]*dlarge[2])
                h11=w1*(rotation*ss)**diagonal_shift*small[:,0]*large[2]-w2*(rotation*t)**diagonal_shift*small[:,1]*large[3]
            else:
                h21=(w1+w2)*small[:,2]*large[1]
                dh21=(w1+w2)*(dsmall[:,2]*large[1]+small[:,2]*dlarge[1])
                h11=-w1*(rotation*ss)**diagonal_shift*small[:,3]*large[1]+w2*(rotation*t)**diagonal_shift*small[:,2]*large[0]
            answer+=rotation**(alpha+beta)*weight*np.array([np.sum(ws*ext*h21),np.sum(ws*ext*h11),np.sum(ws*ext*dh21)])
    # Opposite-SK terms and their signs, integrated on independent contours.
    answer[0]-=lap[1,alpha][2]*lap[-1,beta][1]+lap[-1,alpha][1]*lap[1,beta][2]
    answer[1]+=lap[1,alpha+diagonal_shift][3]*lap[-1,beta][1]-lap[-1,alpha+diagonal_shift][0]*lap[1,beta][2]
    answer[2]-=(dlap[1,alpha][2]*lap[-1,beta][1]+lap[1,alpha][2]*dlap[-1,beta][1]
                +dlap[-1,alpha][1]*lap[1,beta][2]+lap[-1,alpha][1]*dlap[1,beta][2])
    return answer


def seed_values(mu,alpha,beta,diagonal_shift=0):
    nu=mu-.5j
    h21=-sum((1j*a)**j*(1j*b)**k*scalar_seed_folded(a,b,alpha-2+j,beta-2+k,nu)
             for a in (1,-1) for b in (1,-1) for j in (0,1) for k in (0,1))
    # alpha,beta already contain the common positive regulator.
    h11=regulated_diagonal_seed(mu,alpha+diagonal_shift,beta,0,28)
    return np.array([complex(h21),complex(h11)])


if __name__=='__main__':
    lines=[]
    for mu,kind in ((.2,'clock'),(.73,'clock'),(.73,'mixed')):
        delta=.8
        alpha=(1j*mu if kind=='clock' else -1j*mu)+delta
        beta=1j*mu+delta
        shift=int(kind=='clock')
        target=seed_values(mu,alpha,beta,shift)
        coarse=integrate(mu,alpha,beta,48,26.,shift)
        fine=integrate(mu,alpha,beta,72,30.,shift)
        if kind=='clock':
            p1,p2,_,_=regulated_hard_coefficients(mu,delta,0,28)
            direct_p1=fine[0]+2j/(1+2j*mu)*fine[1]
            direct_p2=fine[2]-fine[0]
            perr=max(abs(direct_p1-complex(p1)),abs(direct_p2-complex(p2)))/max(1.,abs(complex(p1)),abs(complex(p2)))
            assert perr<4e-6,(mu,perr)
            lines.append(f'PASS direct regulated first-nonzero hard coefficients P1/P2 mu={mu}: relative error {perr:.3e}')
            print(lines[-1],flush=True)
        fine=fine[:2];coarse=coarse[:2]
        error=np.linalg.norm(fine-target)/max(1.,np.linalg.norm(target))
        stability=np.linalg.norm(fine-coarse)/max(1.,np.linalg.norm(fine))
        assert error<4e-6,(mu,kind,error,fine,target)
        assert stability<2e-5,(mu,kind,stability)
        line=f'PASS direct BD-contour H21/H11 quadrature mu={mu}, {kind}, delta={delta}: relative seed error {error:.3e}; N=48/72, tail=26/30 change {stability:.3e}'
        print(line,flush=True);lines.append(line)
    lines.append('Scope: independently normalized physical-mass time kernels at positive regulator; zero-regulator extrapolation still uses scalar seeds, not this quadrature.')
    Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
