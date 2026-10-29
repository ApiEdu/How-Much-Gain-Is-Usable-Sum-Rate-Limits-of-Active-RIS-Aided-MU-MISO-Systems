import re,collections,numpy as np
S='/tmp/claude-0/-home-claude/001c75b9-0bc8-5d1d-808e-ce0003cea96f/scratchpad/'
def parse(pat,fn,tag):
    out=collections.defaultdict(dict)
    try:
        for ln in open(fn):
            m=re.match(tag+r' (\d+) (\d+) (-?\d+) it (\d+).*?simrate ([\d.]+)',ln)
            if m: out[(int(m.group(1)),int(m.group(2)),int(m.group(3)))][int(m.group(4))]=float(m.group(5))
    except FileNotFoundError: pass
    return out
A=collections.defaultdict(dict);C=collections.defaultdict(dict);U=collections.defaultdict(dict)
for L in (100,200):
    for k,v in parse(None,S+f'convgrad_L{L}.log','CONV').items(): A[k].update(v)
    for k,v in parse(None,S+f'convlin_L{L}.log','CONVLIN').items(): C[k].update(v)
    for k,v in parse(None,S+f'convun_L{L}.log','CONVUN').items(): U[k].update(v)
rows=[];st=collections.defaultdict(list)
for k in sorted(A):
    if not (k in C and k in U): continue
    a,c,u=A[k],C[k],U[k]
    L,d,b=k
    r=(L,b,d,u[50],u[150],c[50],c[150],a[50],a[150])
    rows.append(r)
    st['loss50'].append(c[50]-u[50]); st['loss150'].append(c[150]-u[150])
    st['diff50'].append(a[50]-c[50]); st['diff150'].append(a[150]-c[150])
    st['aware_gain'].append(a[150]-a[50]); st['cap_gain'].append(c[150]-c[50]); st['un_gain'].append(u[150]-u[50])
for r in rows:
    L,b,d,u50,u150,c50,c150,a50,a150=r
    print(f"{L} & ${b:+d}$ & {d} & {u50:.2f} & {u150:.2f} & {c50:.2f} & {c150:.2f} & {a50:.2f} & {a150:.2f} & ${c150-u150:+.2f}$ & ${a50-c50:+.2f}$ & ${a150-c150:+.2f}$\\\\")
for k,v in st.items(): print(k,'min %.2f max %.2f mean %.2f n=%d'%(min(v),max(v),np.mean(v),len(v)))
