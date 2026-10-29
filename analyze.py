import numpy as np, scipy.io as sio, os
bt=[-5,0,5]; names=['FP (linear)','FP + gain cap','Grad-unaware','Grad-linregion','Grad-aware','Grad-aware2']
out={}
for L in (200,100):
    m=sio.loadmat(f'v6/res_main_L{L}.mat'); e=sio.loadmat(f'v6/res_extra_L{L}.mat'); f=sio.loadmat(f'v6/res_fp_L{L}.mat')
    bdb=list(m['bdb'].ravel()); ix=[bdb.index(b) for b in bt]
    R=np.full((14,3,6),np.nan)
    R[:,:,0]=f['R'][:,:,0]; R[:,:,1]=f['R'][:,:,1]
    R[:4,:,2]=e['Run'][:4][:,ix]; R[:4,:,3]=e['Rlr'][:4][:,ix]; R[:4,:,4]=m['Rr'][:4][:,ix,5]
    if os.path.exists(f'v6/res_more_L{L}.mat'):
        g=sio.loadmat(f'v6/res_more_L{L}.mat')['R']; n=int(sio.loadmat(f'v6/res_more_L{L}.mat')['dd'].ravel()[0])
        for s_,k in ((4,0),(3,1),(5,2),(2,3)): R[4:4+n,:,s_]=g[:n,:,k]
    if os.path.exists(f'v6/res_more2_L{L}.mat'):
        g2=sio.loadmat(f'v6/res_more2_L{L}.mat'); n2=int(g2['dd'].ravel()[0])
        for s_,k in ((4,0),(3,1),(5,2),(2,3)): R[11:11+n2,:,s_]=g2['R'][:n2,:,k]
    if os.path.exists(f'v6/res_more3_L{L}.mat'):
        g3=sio.loadmat(f'v6/res_more3_L{L}.mat'); n3=int(g3['dd'].ravel()[0])
        for s_,k in ((4,0),(3,1),(5,2),(2,3)): R[12:12+n3,:,s_]=g3['R'][:n3,:,k]
    out[L]=R
    nd=int(np.sum(~np.isnan(R[:,0,4])))
    print(f'\nL={L} m, draws with all gradient designs: {nd}')
    print(f'{"design":16s}'+''.join(f'  beta={b:+d}dB      ' for b in bt))
    for s,nm in enumerate(names):
        row=''
        for ib in range(3):
            x=R[:,ib,s]; x=x[~np.isnan(x)]
            row+=f'  {x.mean():6.2f}+-{(x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else 0):4.2f} n={len(x):<2d}'
        print(f'{nm:16s}'+row)
    for a,b,lab in ((4,3,'aware - linregion'),(5,3,'aware2 - linregion'),(3,1,'linregion - FP+cap'),(3,0,'linregion - FP')):
        row=''
        for ib in range(3):
            d=R[:,ib,a]-R[:,ib,b]; d=d[~np.isnan(d)]
            row+=f'  {d.mean():+6.2f} (>0:{int((d>0).sum())}/{len(d)})   '
        print(f'{lab:20s}'+row)
np.save('combined_R.npy',np.stack([out[200],out[100]]))
