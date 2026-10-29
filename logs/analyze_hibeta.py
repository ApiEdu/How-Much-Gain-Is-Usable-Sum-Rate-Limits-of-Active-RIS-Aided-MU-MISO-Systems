import sys, collections, math
import numpy as np
rows=collections.defaultdict(list)
for L in (100,200):
    for ln in open(f'hbcf_L{L}.log'):
        t=ln.split()
        if len(t)==6 and t[0]=='HBCF':
            rows[(int(t[1]),int(t[3]))].append((int(t[2]),float(t[4]),float(t[5])))
print('L beta n  unaware  aware  loss(aware-unaware) SE  min max')
for L in (100,200):
    for b in sorted({k[1] for k in rows if k[0]==L}):
        r=rows[(L,b)]
        d=np.array([x[2]-x[1] for x in r]); u=np.array([x[1] for x in r]); a=np.array([x[2] for x in r])
        n=len(d); se=d.std(ddof=1)/math.sqrt(n) if n>1 else float('nan')
        print(f'{L} {b:3d} {n:2d} {u.mean():7.3f} {a.mean():7.3f} {d.mean():+7.3f} {se:6.3f} {d.min():+7.3f} {d.max():+7.3f}')
