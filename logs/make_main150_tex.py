import numpy as np, sys
sys.argv=['x']
import analyze_main150 as A
BT,LS,DR=A.BT,A.LS,A.DR
cols=[(L,b) for L in LS for b in BT]
def cell(des,L,b):
    m,se,n=A.stats([A.val(des,L,d,b) for d in DR]); return m,se,n
rows=[]
for des,nm in A.names:
    cs=[]
    for L,b in cols:
        m,se,n=cell(des,L,b); cs.append(f'${m:.2f}\\pm{se:.2f}$')
    rows.append(nm.replace('J98 / WMMSE (linear)','J98 / WMMSE (linear)')+' & '+' & '.join(cs)+'\\\\')
print('\n'.join(rows))
print('---- paired')
def diffs(a,b,it_a=150,it_b=150):
    out={}
    for L,bb in cols:
        out[(L,bb)]=np.array([A.val(a,L,d,bb,it_a)-A.val(b,L,d,bb,it_b) for d in DR])
    return out
def line(lab,a,b):
    ds=diffs(a,b); cs=[]
    for k in cols:
        x=ds[k]; cs.append(f'${x.mean():+.2f}$ ({(x>0).sum()}/{len(x)})')
    return lab+' & '+' & '.join(cs)+'\\\\'
for lab,a,b in (('Cap $-$ linear model','linreg','unaware'),('Aware $-$ cap','aware','linreg'),('J98 in loop $-$ J98 after','inloop','fpcap'),('J98 in loop $-$ grad. cap','inloop','linreg')):
    print(line(lab,a,b))
print('---- ranges over all individual designs')
def allr(a,b,it_a=150,it_b=150):
    x=np.array([A.val(a,L,d,bb,it_a)-A.val(b,L,d,bb,it_b) for L in LS for bb in BT for d in DR]); x=x[~np.isnan(x)]
    return x
for lab,a,b,ia,ib in (('cap-linear','linreg','unaware',150,150),('cap50-lin50','linreg','unaware',50,50),('aware-cap150','aware','linreg',150,150),('aware-cap50','aware','linreg',50,50),('inloop-cap150','inloop','linreg',150,150),('inloop-posthoc','inloop','fpcap',150,150),('inloop-J98','inloop','fp',150,150)):
    x=allr(a,b,ia,ib); print(f'{lab:16s} n={len(x)} min {x.min():+.2f} max {x.max():+.2f} mean {x.mean():+.2f} pos {(x>0).sum()}')
for des in ('aware','linreg','unaware'):
    x=allr(des,des,150,50); print(f'gain 50->150 {des}: n={len(x)} min {x.min():.2f} max {x.max():.2f} mean {x.mean():.2f}')
x=allr('aware','aware',150,50)
