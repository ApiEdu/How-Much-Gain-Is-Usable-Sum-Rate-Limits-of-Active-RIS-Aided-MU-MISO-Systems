import sys,re,glob,os,numpy as np
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
sys.path.insert(0,S); sys.argv=['x']
import analyze_main150 as A
for fn in sorted(glob.glob(S+'ext150_L*_d*.log')):
    for ln in open(fn):
        m=re.match(r'M150 (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)',ln)
        if m: A.put(m.group(4),int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(5)),float(m.group(6)))
D=range(1,11); LS=(100,200); BT=(-5,0,5)
v=lambda des,L,d,b,it: A.val(des,L,d,b,it)
def ms(x): x=np.asarray(x); return x.mean(), x.std(ddof=1)/np.sqrt(len(x))
if __name__=='__main__':
    print('== table rows (150 it, n=10)')
    for des,nm in (('unaware','lin'),('linreg','cap'),('aware','aware')):
        print(nm,' & '.join('$%.2f\\pm%.2f$'%ms([v(des,L,d,b,150) for d in D]) for L in LS for b in BT))
    for a,b2,nm in (('linreg','unaware','cap-lin'),('aware','linreg','aware-cap')):
        out=[]
        for L in LS:
            for b in BT:
                x=np.array([v(a,L,d,b,150)-v(b2,L,d,b,150) for d in D]); out.append('$%+.2f$ (%d/%d)'%(x.mean(),(x>0).sum(),len(x)))
        print(nm,' & '.join(out))
    print('== ranges over 60 designs')
    def allr(a,b2,ia,ib):
        return np.array([v(a,L,d,b,ia)-v(b2,L,d,b,ib) for L in LS for b in BT for d in D])
    for lab,a,b2,ia,ib in (('cap-lin 150','linreg','unaware',150,150),('cap-lin 50','linreg','unaware',50,50),('aware-cap 150','aware','linreg',150,150),('aware-cap 50','aware','linreg',50,50)):
        x=allr(lab and a,b2,ia,ib); print(lab,'n',len(x),'min %.2f max %.2f mean %.2f pos %d'%(x.min(),x.max(),x.mean(),(x>0).sum()))
    for des in ('aware','linreg','unaware'):
        x=allr(des,des,150,50); print('gain 50->150',des,'min %.2f max %.2f mean %.2f'%(x.min(),x.max(),x.mean()))
    print('== per-condition aware-cap 150 min/max')
    for L in LS:
        for b in BT:
            x=np.array([v('aware',L,d,b,150)-v('linreg',L,d,b,150) for d in D]); print(L,b,'min %+.2f max %+.2f'%(x.min(),x.max()), np.round(x,2))
    print('== per-condition gain 50->150 aware, cap, lin (min..max)')
    for L in LS:
        for b in BT:
            print(L,b,*['%s %.2f..%.2f'%(des,min(v(des,L,d,b,150)-v(des,L,d,b,50) for d in D),max(v(des,L,d,b,150)-v(des,L,d,b,50) for d in D)) for des in ('aware','linreg','unaware')])
    print('== conv table (n=10)')
    for L in LS:
        for b in BT:
            m=lambda des,it: np.mean([v(des,L,d,b,it) for d in D])
            d50=np.mean([v('aware',L,d,b,50)-v('linreg',L,d,b,50) for d in D]); d150=np.mean([v('aware',L,d,b,150)-v('linreg',L,d,b,150) for d in D])
            print('%d & $%+d$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & $%+.2f$ & $%+.2f$\\\\'%(L,b,m('unaware',50),m('unaware',150),m('linreg',50),m('linreg',150),m('aware',50),m('aware',150),d50,d150) if b<0 else '%d & $%s%d$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & $%+.2f$ & $%+.2f$\\\\'%(L,'+' if b>0 else '',b,m('unaware',50),m('unaware',150),m('linreg',50),m('linreg',150),m('aware',50),m('aware',150),d50,d150))
    print('nan check', sum(np.isnan(v(des,L,d,b,it)) for des in ('aware','linreg','unaware') for L in LS for b in BT for d in D for it in (50,150)))
