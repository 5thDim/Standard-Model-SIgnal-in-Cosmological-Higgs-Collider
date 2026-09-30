"""Symmetry and engineering-dimension checks for the final result evaluator.
These are consistency checks of the assembled formula, not independent loop
integrations. Cache only mass-dependent scalar functions across geometries.
"""
from functools import lru_cache
from pathlib import Path
import mpmath as mp
import top_box_signal as signal
mp.mp.dps=35
signal.clock_angular_coefficients=lru_cache(None)(signal.clock_angular_coefficients)
signal.mixed_hard_seed=lru_cache(None)(signal.mixed_hard_seed)
L=mp.matrix([mp.mpf('.8'),0,mp.mpf('.6')])
R=mp.matrix([mp.mpf('.3'),mp.mpf('1.2'),mp.mpf('-.4')])
S=mp.matrix([0,0,mp.mpf('.02')])
k=[L+S/2,-L+S/2,R-S/2,-R-S/2]
base=signal.collapsed_signal(k,'.73')
lines=[]
for name,momenta,H,scale in [
    ('left Bose exchange',[k[1],k[0],k[2],k[3]],1,1),
    ('right Bose exchange',[k[0],k[1],k[3],k[2]],1,1),
    ('hard-pair exchange',[k[2],k[3],k[0],k[1]],1,1),
    ('parity',[-v for v in k],1,1),
    ('dimension -5 including H',[2*v for v in k],2,mp.mpf(2)**-5),
]:
    result=signal.collapsed_signal(momenta,'.73',H=H)
    residual=max(abs(result[field]/(base[field]*scale)-1) for field in ('clock','mixed'))
    assert residual<mp.mpf('1e-28'),(name,residual)
    lines.append(f'PASS {name}: relative residual {float(residual):.3e}')
lines.append('Scope: symmetries and dimensions of the assembled leading formula only.')
Path(__file__).with_suffix('.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
