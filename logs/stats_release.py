import re, collections
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
DJ={}; SC={}; PST={}
for L in (100,200):
    for ln in open(S+'warm_L%d.log'%L):
        m=re.match(r'DJ (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)',ln)
        if m: DJ[(int(m[1]),int(m[3]),int(m[2]),m[4],int(m[5]))]=tuple(float(m[i]) for i in (6,7,8))
        m=re.match(r'WARM (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)',ln)
        if m: SC[(int(m[1]),int(m[3]),int(m[2]),m[4],int(m[5]))]=float(m[6])
        m=re.match(r'PST (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)',ln)
        if m: PST[(int(m[1]),int(m[3]),int(m[2]),m[4],int(m[5]))]=tuple(float(m[i]) for i in (6,7,8))
keys=[(L,b,d) for L in (100,200) for b in (-5,0,5) for d in (1,2,3,4)]
print('designs',len(keys), len(DJ), len(SC), len(PST))
print('kept50 J > start J:', sum(DJ[k+('capcont',50)][0]>DJ[k+('cap',50)][0] for k in keys), '/24')
print('kept100 J > start J:', sum(DJ[k+('capcont',100)][0]>DJ[k+('cap',50)][0] for k in keys))
print('released100 J > start J:', sum(DJ[k+('awarewarm',100)][0]>DJ[k+('cap',50)][0] for k in keys))
print('released50 J == start J:', sum(DJ[k+('awarewarm',50)][0]==DJ[k+('cap',50)][0] for k in keys), [k for k in keys if DJ[k+('awarewarm',50)][0]==DJ[k+('cap',50)][0]])
print('released100 J == start J:', sum(DJ[k+('awarewarm',100)][0]==DJ[k+('cap',50)][0] for k in keys))
d=[SC[k+('awarewarm',100)]-SC[k+('capcont',100)] for k in keys]
print('released-kept symbol: min %.2f max %.2f mean %.3f pos %d'%(min(d),max(d),sum(d)/24,sum(x>0 for x in d)))
d2=[SC[k+('awarewarm',100)]-SC[k+('cap',50)] for k in keys]
print('released100-start: min %.2f max %.2f'%(min(d2),max(d2)), 'kept100-start', min(SC[k+('capcont',100)]-SC[k+('cap',50)] for k in keys), max(SC[k+('capcont',100)]-SC[k+('cap',50)] for k in keys))
for L in (100,200):
    sub=[k for k in keys if k[0]==L]
    dd=[SC[k+('awarewarm',100)]-SC[k+('capcont',100)] for k in sub]
    print(L,'min %.2f max %.2f pos %d/12'%(min(dd),max(dd),sum(x>0 for x in dd)))
for tag,it in (('capcont',100),('awarewarm',100)):
    v=[PST[k+(tag,it)] for k in keys]
    print(tag,'max p (max over designs) %.1f, frac>=0.99pmax max %.3f mean %.4f'%(max(x[1] for x in v),max(x[2] for x in v),sum(x[2] for x in v)/24),
          'designs with any element >=0.99pmax:', sum(x[2]>0 for x in v))
