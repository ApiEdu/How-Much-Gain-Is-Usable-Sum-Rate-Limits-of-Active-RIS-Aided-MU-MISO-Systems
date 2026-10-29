# Regenerates fig_sumrate_vs_beta.pdf from the numbers printed in Tables VI and VII of paper.tex (means and standard errors over 14 draws).
import re,sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
t=open(S+'paper/paper.tex',encoding='utf-8').read()
def table(label):
    i=t.index(r'\label{'+label+'}'); j=t.index(r'\end{tabular}',i)
    rows=[r for r in t[i:j].split('\n') if r.endswith(r'\\') and r'\pm' in r]
    out=[]
    for r in rows:
        name=r.split('&')[0].strip()
        vals=[(float(a),float(b)) for a,b in re.findall(r'\$([\d.]+)\\pm([\d.]+)\$',r)]
        out.append((name,vals))
    return out
T={100:table('tab:L100'),200:table('tab:L200')}
assert all(len(v)==6 for v in T.values()) and all(len(x[1])==3 for v in T.values() for x in v)
labels=['WMMSE (linear)','WMMSE + gain cap','Gradient, linear model','Gradient, linear-region cap','Gradient, compression-aware','Gradient, compression-aware (start: cap)']
style=[('#666666','s','--'),('#666666','D','-'),('#D55E00','v','--'),('#009E73','o','-'),('#8172B3','^','-'),('#E7298A','*','-')]
plt.rcParams.update({'pdf.fonttype':42,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(6.1,3.3))
beta=[-5,0,5]
for ax,L in zip(axs,(100,200)):
    for (name,vals),lab,(c,m,ls) in zip(T[L],labels,style):
        ax.errorbar(beta,[v[0] for v in vals],yerr=[v[1] for v in vals],color=c,marker=m,ls=ls,lw=1.6,ms=5,capsize=2,label=lab)
    ax.set_title(f'$L={L}$ m'); ax.set_xlabel(r'Saturation ratio $\beta$ (dB)'); ax.set_xticks(beta); ax.grid(alpha=.3)
axs[0].set_ylabel('Sum-rate (bit/s/Hz)')
h,l=axs[0].get_legend_handles_labels()
fig.legend(h,l,loc='lower center',ncol=2,frameon=False,fontsize=7.5)
fig.tight_layout(rect=(0,0.2,1,1))
fig.savefig(S+'mdpi_build/fig_sumrate_vs_beta.pdf'); fig.savefig(S+'mdpi_build/fig_sumrate_vs_beta.png',dpi=600)
print('ok')
