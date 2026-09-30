"""Evaluate leading collapsed top-box terms in Top Quark Signal.tex.

The result is for two identical Majorana copies and Nc colors, with massless
external Higgs contractions. It is not Aoki's bubble amplitude. Contact terms
and subleading soft terms are excluded. Momentum units must match H.

Example: python3 checks/top_box_signal.py --mu .73 --soft .02
"""
import argparse
import json
import mpmath as mp
from weyl_scalar_seed import clock_angular_coefficients, mixed_hard_seed


def collapsed_signal(momenta, mu, gt=1, H=1, Nc=3, dps=35):
    """Return clock and mixed contributions, with delta function removed.

    k1+k2 defines the collapsed channel. Caller should check soft_ratio << 1.
    Overall engineering dimension is -5. The integral retains complete hard
    SK sums; analytic local terms have been removed by momentum continuation.
    """
    with mp.workdps(dps):
        ks=[mp.matrix([mp.mpf(v) for v in k]) for k in momenta]
        if len(ks)!=4 or any(len(k)!=3 for k in ks):
            raise ValueError('Expected four three-component momenta')
        norms=[mp.norm(k) for k in ks]
        if min(norms)<=0: raise ValueError('External momenta must be nonzero')
        if mp.norm(sum(ks,mp.matrix([0,0,0])))>mp.mpf('1e-12')*max(norms):
            raise ValueError('External momenta must conserve momentum')
        soft=ks[0]+ks[1]; s=mp.norm(soft)
        L=(ks[0]-ks[1])/2;R=(ks[2]-ks[3])/2
        kL,kR=mp.norm(L),mp.norm(R)
        if min(s,kL,kR)<=0: raise ValueError('Use nonzero soft and hard momenta')
        mu,gt,H=map(mp.mpf,(mu,gt,H))
        if mu<=0: raise ValueError('Separated branches require mu>0; massless limit must be combined first')
        n,m,z=L/kL,R/kR,soft/s
        c1,c3=mp.fdot(n,z),mp.fdot(m,z)
        transverse=mp.fdot(n,m)-c1*c3
        # Polynomial angular invariants avoid ill-defined azimuths at a pole.
        coefficients=clock_angular_coefficients(mu,dps)
        A1,A2,A3,B1,D1=[coefficients[k] for k in ('A1','A2','A3','B1','D1')]
        angular=(A1+A2*(c1**2+c3**2)+A3*c1**2*c3**2
                 +4*B1*c1*c3*transverse
                 +D1*(2*transverse**2-(1-c1**2)*(1-c3**2)))
        prefactor=Nc*gt**4*H**4/(16*mp.fprod(k**3 for k in norms))
        ratio=s**2/(4*(norms[0]+norms[1])*(norms[2]+norms[3]))
        clock=-prefactor*s**5/(kL*kR)*2*mp.re(ratio**(2j*mu)*angular)
        mixed_h=mixed_hard_seed(mu,dps)
        mixed=(prefactor*s**3*2*mu*(1+mu**2)/(3*mp.pi*mp.sinh(2*mp.pi*mu))
               *mp.re(mixed_h)**2)
        return dict(clock=+clock,mixed=+mixed,total=+(clock+mixed),
                    soft_ratio=+max(s/kL,s/kR),angular=angular,
                    mixed_seed_imaginary_residual=+abs(mp.im(mixed_h)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mu',default='.73')
    parser.add_argument('--soft',default='.02')
    args=parser.parse_args()
    with mp.workdps(35):
        s=mp.mpf(args.soft);L=mp.matrix([mp.mpf('.8'),0,mp.mpf('.6')])
        R=mp.matrix([mp.mpf('.3'),mp.mpf('1.2'),mp.mpf('-.4')])
        S=mp.matrix([0,0,s])
        result=collapsed_signal([L+S/2,-L+S/2,R-S/2,-R-S/2],args.mu)
    print(json.dumps({k:mp.nstr(v,18) for k,v in result.items()},indent=2))
