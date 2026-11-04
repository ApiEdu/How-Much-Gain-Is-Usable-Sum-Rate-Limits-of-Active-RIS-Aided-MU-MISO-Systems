"""Equal-budget (150-iteration) gradient comparison over up to 14 channel draws.
Draws 1-4 come from main150_L*.log (+ earlier conv logs, via analyze_main150), draws 5-14 from ext150_L*_d*.log.
Usage: python3 -I analyze14.py"""
import re, glob, sys, numpy as np
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
sys.path.insert(0,S); sys.argv=['x']
import analyze_main150 as A
for fn in sorted(glob.glob(S+'ext150_L*_d*.log')):
    for ln in open(fn):
        m=re.match(r'M150 (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)',ln)
        if m: A.put(m.group(4),int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(5)),float(m.group(6)))
BT,LS=A.BT,A.LS
def have(L,b,dr):  # draws with all three designs at 150
    return [d for d in dr if all(not np.isnan(A.val(x,L,d,b,150)) for x in ('unaware','linreg','aware'))]
def table(dr_all=range(1,15)):
    out={}
    for L in LS:
        for b in BT:
            ds=have(L,b,dr_all)
            out[(L,b)]=ds
    return out
if __name__=='__main__':
    t=table()
    print('draws available per condition:',{k:len(v) for k,v in t.items()})
    for des in ('unaware','linreg','aware'):
        row=[]
        for L in LS:
            for b in BT:
                m,se,n=A.stats([A.val(des,L,d,b,150) for d in t[(L,b)]]); row.append(f'{m:6.2f}+-{se:4.2f}(n{n})')
        print(f'{des:8s}',' | '.join(row))
    for a,bn,lab in (('linreg','unaware','cap - linear'),('aware','linreg','aware - cap')):
        row=[]
        for L in LS:
            for b in BT:
                x=np.array([A.val(a,L,d,b,150)-A.val(bn,L,d,b,150) for d in t[(L,b)]])
                row.append(f'{x.mean():+6.2f}+-{x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else float("nan"):4.2f}({(x>0).sum()}/{len(x)})')
        print(f'{lab:13s}',' | '.join(row))
