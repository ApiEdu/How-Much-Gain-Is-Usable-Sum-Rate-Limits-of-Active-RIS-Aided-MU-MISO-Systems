import numpy as np, collections, glob, sys
rows=[]
for f in sorted(glob.glob('/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/valid_L*.log')):
    for ln in open(f):
        if ln.startswith('VAL '):
            p=ln.split(); rows.append((p[1],int(p[2]),int(p[3]),int(p[4]),int(p[5]),float(p[6]),float(p[7]),float(p[8])))
print(len(rows),'rows')
def rate(x): return np.log2(1+10**(np.asarray(x)/10))
# per-design (group,L,d,beta): rate from users
des=collections.defaultdict(list)
for g,L,d,b,k,m,i,s in rows: des[(g,L,d,b)].append((m,i,s))
out=[]
for grp in ('CF-aware','CF-unaware','Grad-aware'):
  for b in (-5,0,5):
    sel=[(m,i,s) for g,L,d,bb,k,m,i,s in rows if g==grp and bb==b]
    if not sel: continue
    a=np.array(sel); em=a[:,0]-a[:,2]; ei=a[:,1]-a[:,2]
    ds=[v for (g,L,d,bb),v in des.items() if g==grp and bb==b]
    rm=[np.sum(rate([x[0] for x in v]))-np.sum(rate([x[2] for x in v])) for v in ds]
    ri=[np.sum(rate([x[1] for x in v]))-np.sum(rate([x[2] for x in v])) for v in ds]
    print(f'{grp:11s} b={b:+d} n_users={len(a):3d} ndes={len(ds):2d} | corr: mean {em.mean():+.2f} mean|e| {np.abs(em).mean():.2f} max|e| {np.abs(em).max():.2f} dRate {np.mean(rm):+.2f} | indep: mean {ei.mean():+.2f} mean|e| {np.abs(ei).mean():.2f} max|e| {np.abs(ei).max():.2f} dRate {np.mean(ri):+.2f}')
    out.append((grp,b,len(ds),np.abs(em).mean(),np.abs(em).max(),np.mean(rm),np.abs(ei).mean(),np.abs(ei).max(),np.mean(ri)))
a=np.array([(m-s,i-s) for g,L,d,b,k,m,i,s in rows]); print('ALL corr mean|e| %.2f max %.2f | indep mean|e| %.2f max %.2f'%(np.abs(a[:,0]).mean(),np.abs(a[:,0]).max(),np.abs(a[:,1]).mean(),np.abs(a[:,1]).max()))
