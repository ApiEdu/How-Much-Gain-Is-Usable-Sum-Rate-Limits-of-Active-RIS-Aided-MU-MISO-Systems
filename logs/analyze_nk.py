import re,collections,numpy as np,glob
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
D=collections.defaultdict(dict)
for f in sorted(glob.glob(S+'convnk_N128_K*_L*.log')):
    for ln in open(f):
        m=re.match(r'NKEQ N=128 K=(\d) L=(\d+) d=(\d+) b=(-?\d+) (\w+) it (\d+) simrate ([\d.]+)',ln)
        if m:
            K,L,d,b=map(int,m.groups()[:4]); D[(K,L,d,b)][(m.group(5),int(m.group(6)))]=float(m.group(7))
all50=[];all150=[];rows=[]
for k in sorted(D):
    v=D[k]
    if len(v)<6: continue
    K,L,d,b=k; d50=v[('aware',50)]-v[('linreg',50)]; d150=v[('aware',150)]-v[('linreg',150)]
    all50.append(d50); all150.append(d150)
    print(k,'aware %.2f/%.2f cap %.2f/%.2f  d50 %+.2f d150 %+.2f'%(v[('aware',50)],v[('aware',150)],v[('linreg',50)],v[('linreg',150)],d50,d150))
    rows.append((K,L,b,d,d50,d150))
print('n',len(all150),'d50 %.2f..%.2f  d150 %.2f..%.2f mean150 %.2f'%(min(all50),max(all50),min(all150),max(all150),np.mean(all150)))
# per-config means
for K in (2,3):
    for L in (100,200):
        for b in (-5,5):
            x=[r[5] for r in rows if r[0]==K and r[1]==L and r[2]==b]
            y=[r[4] for r in rows if r[0]==K and r[1]==L and r[2]==b]
            if x: print(f'K={K} L={L} b={b:+d}: d50 mean {np.mean(y):+.2f}  d150 mean {np.mean(x):+.2f} (n={len(x)})')
