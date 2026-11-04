import re, sys, glob, json
sys.argv=['x']
import analyze_main150 as A
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
R={}; P={}
for f in glob.glob(S+'rho_L*.log')+glob.glob(S+'pmax*_L*.log'):
    for ln in open(f):
        m=re.match(r'RHO (\d+) (\d+) (-?\d+) ([\d.]+) (\d+) (\d+) ([\d.]+)',ln)
        if m: R[(int(m[1]),float(m[4]),int(m[5]),int(m[2]))]=float(m[7])
        m=re.match(r'PST (\d+) (\d+) (-?\d+) rho([\d.]+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)',ln)
        if m: P[(int(m[1]),float(m[4]),int(m[2]),'f'+f.split('/')[-1].split('_')[0])]=(float(m[6]),float(m[7]),float(m[8]))
for L in (100,200):
    for d in (1,2,3): R[(L,0.1,100,d)]=A.val('linreg',L,d,0,150)
aw={(L,d):A.val('aware',L,d,0,150) for L in (100,200) for d in (1,2,3)}
out={}
def r(L,rho,pm,d): return R[(L,rho,pm,d)]
for L in (100,200):
    D=(1,2,3)
    d003=[r(L,0.03,100,d)-r(L,0.1,100,d) for d in D]
    d03=[r(L,0.3,100,d)-r(L,0.1,100,d) for d in D]
    d1=[r(L,1.0,100,d)-r(L,0.1,100,d) for d in D]
    p30=[r(L,0.1,30,d)-r(L,0.1,100,d) for d in D]
    p300=[r(L,0.1,300,d)-r(L,0.1,100,d) for d in D]
    best=[max((0.03,0.1,0.3,1.0),key=lambda x:r(L,x,100,d)) for d in D]
    am=[aw[(L,d)]-r(L,0.1,100,d) for d in D]
    out[L]=dict(d003=d003,d03=d03,d1=d1,p30=p30,p300=p300,best=best,aw_minus_cap=am,
                aw_minus_03=[aw[(L,d)]-r(L,0.3,100,d) for d in D],aw_minus_1=[aw[(L,d)]-r(L,1.0,100,d) for d in D])
    print(L,{k:[round(x,2) if isinstance(x,float) else x for x in v] for k,v in out[L].items()})
for k,v in sorted(P.items(), key=lambda kv: str(kv[0])): print(k,v)
json.dump({str(k):v for k,v in out.items()},open(S+'rho_numbers.json','w'),indent=1)
