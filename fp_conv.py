# Does the J98-style WMMSE/FP reference keep improving with more iterations? Run on 3 draws per distance, log rate at checkpoints.
import sys; sys.argv=['x']; sys.path.insert(0,'.')
import numpy as np, scipy.io as sio, fp_j98 as F
for L in (200,100):
    for d in (1,2,3):
        m=sio.loadmat(f'ch/ch_L{L}_d{d}.mat'); (r,psi,W),hist=F.solve(m['h'],m['f'],m['G'],iters=6000)
        cp=[500,1500,3000,6000]
        print(L,d,' '.join(f'it{c}:{max(hist[:c+1]):.2f}' for c in cp),f'Pris={F.ris_power(m["G"],psi,W)/F.Pris:.2f}',flush=True)
