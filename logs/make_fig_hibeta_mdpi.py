# Loss of the compression-blind design versus the saturation ratio beta, up to +30 dB.
#  closed-form: aware minus unaware (14 draws, symbol-level), logs hbcf_L*.log
#  gradient   : capped minus linear-model design after 150 iterations (draws 1-4): beta -5/0/+5 from main150 logs, beta 10/20 from hbg logs
import re, sys, collections, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'.')
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
import analyze_main150 as A
def se(x):
    x=np.asarray(x,float); return x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else np.nan
cf={}
for L in (100,200):
    rows=collections.defaultdict(dict)
    for ln in open(S+f'hbcf_L{L}.log'):
        t=ln.split()
        if len(t)==6 and t[0]=='HBCF': rows[int(t[3])][int(t[2])]=float(t[5])-float(t[4])
    cf[L]=rows
hbg=collections.defaultdict(dict)
for L in (100,200):
    try:
        for ln in open(S+f'hbg_L{L}.log'):
            t=ln.split()
            if len(t)==7 and t[0]=='HBG': hbg[(L,t[4],int(t[1]) if False else int(t[2]),int(t[3]),int(t[5]))]=float(t[6])
    except FileNotFoundError: pass
def grad_loss(L,b):
    out=[]
    for d in (1,2,3,4):
        if b in (-5,0,5): out.append(A.val('linreg',L,d,b)-A.val('unaware',L,d,b))
        else:
            u=hbg.get((L,'unaware',d,b,150)); c=hbg.get((L,'linreg',d,b,150))
            if u is not None and c is not None: out.append(c-u)
    return np.array(out)
plt.rcParams.update({'pdf.fonttype':42,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(6.1,3.3))
res={}
for ax,L in zip(axs,(100,200)):
    bs=sorted(cf[L]); m=[np.mean(list(cf[L][b].values())) for b in bs]; e=[se(list(cf[L][b].values())) for b in bs]
    ax.errorbar(bs,m,yerr=e,color='#8172B3',marker='^',ls='-',lw=1.6,ms=5,capsize=2,label='Closed form: aware $-$ unaware (14 draws)')
    gb=[-5,0,5,10,20]; gm=[];ge=[];gx=[];gn=[]
    for b in gb:
        x=grad_loss(L,b)
        if len(x)>0: gx.append(b); gm.append(x.mean()); ge.append(se(x)); gn.append(len(x))
    ax.errorbar(gx,gm,yerr=ge,color='#009E73',marker='o',ls='--',lw=1.6,ms=5,capsize=2,label='Gradient: cap $-$ linear model, 150 it. (draws 1–4)')
    res[L]=(gx,gm,ge,gn)
    ax.axhline(0,color='k',lw=.6); ax.set_title(f'$L={L}$ m'); ax.set_xlabel(r'Saturation ratio $\beta$ (dB)'); ax.set_xticks([-5,0,5,10,15,20,25,30]); ax.grid(alpha=.3)
axs[0].set_ylabel('Sum-rate gain (bit/s/Hz)')
h,l=axs[0].get_legend_handles_labels()
fig.legend(h,l,loc='lower center',ncol=1,frameon=False,fontsize=7.5)
fig.tight_layout(rect=(0,0.18,1,1))
fig.savefig(S+'mdpi_build/fig_hibeta.pdf'); fig.savefig(S+'mdpi_build/fig_hibeta.png',dpi=600)
for L in res: print(L,[ (b,round(m,2),round(s,2),n) for b,m,s,n in zip(*res[L])])
