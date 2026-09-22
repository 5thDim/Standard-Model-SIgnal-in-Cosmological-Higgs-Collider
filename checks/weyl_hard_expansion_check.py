"""Check the hard expansion on a finite time interval with symmetric quadrature.
This validates momentum derivatives and endpoint signs, not the BD hard-seed value.
"""
import numpy as np
import mpmath as mp
mp.mp.dps=25
mu=.73;d=1j/(1+2j*mu)
sigma=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
def slash(v):return np.einsum('i,ijk->jk',v,sigma)
z,w=np.polynomial.legendre.leggauss(32);ts=.2+(z+1)*1.4;ws=w*1.4
x=ts[:,None];y=ts[None,:];weights=ws[:,None]*ws[None,:]*(x*y)**(1j*mu-1)
th=(x<y).astype(float)+.5*(x==y)
def radial(p,e1=1.,e2=1.):
 f=np.array([complex(mp.sqrt(mp.pi*p*t)/2*mp.exp(mp.pi*mu/2)*mp.hankel1(.5-1j*mu,p*t)) for t in ts]);g=np.array([complex(mp.sqrt(mp.pi*p*t)/2*mp.exp(-mp.pi*mu/2)*mp.hankel1(.5+1j*mu,p*t)) for t in ts])
 gt=[np.outer(f,f.conj()),np.outer(f,g.conj()),np.outer(g,f.conj()),np.outer(g,g.conj())]
 lt=[-np.outer(g.conj(),g),np.outer(g.conj(),f),np.outer(f.conj(),g),-np.outer(f.conj(),f)]
 out=[np.zeros_like(weights) for _ in range(4)]
 for a in [-1,1]:
  for b in [-1,1]:
   ex=(1+1j*a*e1*x)*np.exp(-1j*a*e1*x)*(1+1j*b*e2*y)*np.exp(-1j*b*e2*y)
   for k in range(4):
    r=gt[k] if (a,b)==(-1,1) else lt[k] if (a,b)==(1,-1) else (th*gt[k]+(1-th)*lt[k]) if a==1 else ((1-th)*gt[k]+th*lt[k])
    out[k]+=a*b*weights*ex*r
 return out
r=radial(1);f=r[2].sum();D=(x*r[0]).sum();step=1e-4
fp=(radial(1+step)[2].sum()-radial(1-step)[2].sum())/(2*step)
a=f+2*d*D;b=fp-f
n=np.array([np.sin(.7),0,np.cos(.7)])
errs=[]
for scale in [1e-3,5e-4]:
 ks=scale*np.array([0,0,1.]);q=scale*np.array([.3,-.2,.1]);K=q+ks;ell=q+ks/2
 got=np.zeros((2,2),complex)
 for direction in [1,-1]:
  p=direction*n+ell;pabs=np.linalg.norm(p)
  e1=np.linalg.norm(direction*n+ks/2);e2=np.linalg.norm(-direction*n+ks/2)
  rr=radial(pabs,e1,e2)
  got+=-rr[2].sum()*slash(p/pabs)-d*(x*rr[0]).sum()*slash(q)+d*(y*rr[3]).sum()*slash(K)+d*d*(x*y*rr[1]).sum()*slash(q)@slash(p/pabs)@slash(K)
 expected=-2*(a*slash(ell)+b*(n@ell)*slash(n))
 err=np.linalg.norm(got-expected)/np.linalg.norm(expected)
 errs.append(err)
 print('scale',scale,'relative projected hard error',err)
assert errs[1]<errs[0]/3
assert errs[1]<1e-6
