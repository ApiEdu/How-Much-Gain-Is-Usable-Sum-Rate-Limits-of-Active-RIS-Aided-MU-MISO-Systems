import numpy as np, collections, os, sys
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
def se(x): x=np.asarray(x,float); return x.std(ddof=1)/np.sqrt(len(x))
ref=collections.defaultdict(dict); csi=collections.defaultdict(dict)
for L in (100,200):
    for ln in open(S+f'hbcf_L{L}.log'):
        t=ln.split()
        if len(t)==6 and t[0]=='HBCF': ref[(L,int(t[3]))][int(t[2])]=(float(t[4]),float(t[5]))
    for ln in open(S+f'csi_L{L}.log'):
        t=ln.split()
        if len(t)==7 and t[0]=='CSI': csi[(L,int(t[4]),int(t[3]))][int(t[2])]=(float(t[5]),float(t[6]))
rows=[]; out=[]
for L in (100,200):
    for b in (-5,0,5):
        cells_a=[];cells_g=[]
        ds=sorted(ref[(L,b)])
        pa=np.array([ref[(L,b)][d][1] for d in ds]); pg=np.array([ref[(L,b)][d][1]-ref[(L,b)][d][0] for d in ds])
        cells_a.append(f'{pa.mean():.2f}'); cells_g.append(f'${pg.mean():+.2f}\\pm{se(pg):.2f}$')
        for nm in (-20,-10):
            r=csi[(L,b,nm)]; ds=sorted(r)
            a=np.array([r[d][1] for d in ds]); g=np.array([r[d][1]-r[d][0] for d in ds])
            cells_a.append(f'{a.mean():.2f}'); cells_g.append(f'${g.mean():+.2f}\\pm{se(g):.2f}$')
        out.append(f'{L} & ${b:+d}$ & '+' & '.join(cells_a)+' & '+' & '.join(cells_g)+'\\\\')
open(S+'csi_cf_rows.tex','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
wins=0;tot=0
for k,r in csi.items():
    for d,(u,a) in r.items(): tot+=1; wins+=a>u
print('aware>unaware',wins,tot)
