import re,collections,numpy as np,scipy.io as sio,sys
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
BT=(-5,0,5); LS=(100,200); DR=(1,2,3,4)
G=collections.defaultdict(dict)       # (design,L,d,b) -> {it: rate}
def put(des,L,d,b,it,r): G[(des,L,d,b)][it]=r
for L in LS:
    try:
        for ln in open(S+f'main150_L{L}.log'):
            m=re.match(r'M150 (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)',ln)
            if m: put(m.group(4),int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(5)),float(m.group(6)))
    except FileNotFoundError: pass
    # earlier equal-budget runs (draws 1-2, beta -5/+5)
    for fn,tag,des in ((f'convgrad_L{L}.log','CONV','aware'),(f'convlin_L{L}.log','CONVLIN','linreg')):
        for ln in open(S+fn):
            m=re.match(tag+r' (\d+) (\d+) (-?\d+) it (\d+).*?simrate ([\d.]+)',ln)
            if m and int(m.group(4)) in (50,150): put(des,int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(4)),float(m.group(5)))
fp={L:sio.loadmat(S+f'sim/v6/res_fp6k_L{L}.mat')['R'] for L in LS}
inl={}
try:
    for ln in open(S+'fpcap_scores.txt'):
        p=ln.split(); inl[(int(p[1]),int(p[2]),int(p[3]))]=float(p[5])
except FileNotFoundError: pass
def val(des,L,d,b,it=150):
    if des=='fp': return fp[L][d-1,BT.index(b),0]
    if des=='fpcap': return fp[L][d-1,BT.index(b),1]
    if des=='inloop': return inl.get((L,d,b),np.nan)
    return G[(des,L,d,b)].get(it,np.nan)
names=[('fp','J98 / WMMSE (linear)'),('fpcap','J98 + cap, after'),('inloop','J98 + cap, in loop'),('unaware','Grad., linear model'),('linreg','Grad., linear-region cap'),('aware','Grad., compression-aware')]
def stats(x):
    x=np.asarray(x,float); x=x[~np.isnan(x)]; return (x.mean(), x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else np.nan, len(x))
if __name__=='__main__':
    for des,nm in names:
        row=[]
        for L in LS:
            for b in BT:
                m,se,n=stats([val(des,L,d,b) for d in DR]); row.append(f'{m:6.2f}+-{se:4.2f}(n{n})')
        print(f'{nm:26s}',' | '.join(row))
    print('--- paired differences (mean, wins/n) ---')
    for a,bn,lab in (('aware','linreg','aware - cap'),('linreg','unaware','cap - linear'),('inloop','linreg','in-loop WMMSE - grad cap'),('inloop','fpcap','in-loop - post-hoc'),('aware','inloop','aware - in-loop WMMSE')):
        row=[]
        for L in LS:
            for b in BT:
                x=np.array([val(a,L,d,b)-val(bn,L,d,b) for d in DR]); x=x[~np.isnan(x)]
                row.append(f'{x.mean():+6.2f}({(x>0).sum()}/{len(x)})' if len(x) else '  n/a')
        print(f'{lab:26s}',' | '.join(row))
