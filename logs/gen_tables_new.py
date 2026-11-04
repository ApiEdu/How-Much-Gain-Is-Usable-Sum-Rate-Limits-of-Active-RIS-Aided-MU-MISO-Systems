#!/usr/bin/env python3
"""Generate LaTeX table bodies for the release-the-cap test and the rho_max / p_max sweeps straight from the logs,
and dump all numbers used in the text to new_numbers.json (so the prose can be checked against it)."""
import re, math, json, collections, sys
sys.argv = ['x']
import analyze_main150 as A      # loader of the 150-iteration reference runs (draws 1-4 from several logs)
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
def mean_se(x):
    n = len(x); m = sum(x) / n
    if n < 2: return m, float('nan')
    return m, math.sqrt(sum((a - m) ** 2 for a in x) / (n - 1) / n)

W = collections.defaultdict(dict)
for L in (100, 200):
    for ln in open(S + 'warm_L%d.log' % L):
        m = re.match(r'WARM (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
        if m: W[(int(m[1]), int(m[3]))].setdefault((m[4], int(m[5])), {})[int(m[2])] = float(m[6])
M = collections.defaultdict(dict)
for L in (100, 200):
    for ln in open(S + 'main150_L%d.log' % L):
        m = re.match(r'M150 (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
        if m: M[(int(m[1]), int(m[3]), m[4], int(m[5]))][int(m[2])] = float(m[6])

DJ = {}
for L in (100, 200):
    for ln in open(S + 'warm_L%d.log' % L):
        m = re.match(r'DJ (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)', ln)
        if m: DJ[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = (float(m[6]), float(m[7]), float(m[8]))
conds = [(100, -5), (100, 0), (100, 5), (200, -5), (200, 0), (200, 5)]
num = {}
cells = collections.defaultdict(list)
for c in conds:
    d = W[c]
    ds = sorted(set(d.get(('cap', 50), {})) & set(d.get(('capcont', 100), {})) & set(d.get(('awarewarm', 100), {})))
    num[str(c)] = {'draws': ds}
    if not ds: continue
    st = [d[('cap', 50)][i] for i in ds]
    ck = [d[('capcont', 100)][i] for i in ds]
    aw = [d[('awarewarm', 100)][i] for i in ds]
    df = [a - b for a, b in zip(aw, ck)]
    cold_c = [A.val('linreg', c[0], i, c[1], 150) for i in ds]
    cold_a = [A.val('aware', c[0], i, c[1], 150) for i in ds]
    dc = [a - b for a, b in zip(aw, cold_c)]
    num[str(c)].update(start=mean_se(st), kept=mean_se(ck), released=mean_se(aw), diff=mean_se(df), npos=sum(1 for x in df if x > 0),
                       cold_cap=mean_se(cold_c), cold_aware=mean_se(cold_a), rel_minus_coldcap=mean_se(dc), per_draw_diff=df,
                       diff_cold_aware_minus_cold_cap=mean_se([a - b for a, b in zip(cold_a, cold_c)]))
    L_, b_ = c
    dJ = [DJ[(L_, b_, i, 'awarewarm', 100)][0] - DJ[(L_, b_, i, 'capcont', 100)][0] for i in ds]
    dR = [DJ[(L_, b_, i, 'awarewarm', 100)][1] - DJ[(L_, b_, i, 'capcont', 100)][1] for i in ds]
    gk = [DJ[(L_, b_, i, 'capcont', 100)][0] - DJ[(L_, b_, i, 'capcont', 100)][1] for i in ds]
    ga = [DJ[(L_, b_, i, 'awarewarm', 100)][0] - DJ[(L_, b_, i, 'awarewarm', 100)][1] for i in ds]
    num[str(c)].update(dJ=mean_se(dJ), dR=mean_se(dR), gap_kept=mean_se(gk), gap_rel=mean_se(ga),
                       posJ=sum(1 for x in dJ if x > 0), posR=sum(1 for x in dR if x > 0),
                       unch50=sum(1 for i in ds if DJ[(L_, b_, i, 'awarewarm', 50)][0] == DJ[(L_, b_, i, 'cap', 50)][0]),
                       dJ_list=dJ, dR_list=dR)
def fmt(ms): return '$%.2f\\pm%.2f$' % ms
def fmt_d(ms, n, npos): return '$%+.2f$ (%d/%d)' % (ms[0], npos, n)
rows = []
def row(label, key, f=fmt):
    cs = []
    for c in conds:
        v = num[str(c)]
        cs.append(f(v[key]) if key in v else '--')
    rows.append(label + ' & ' + ' & '.join(cs) + '\\\\')
row('Capped design, 50 iterations (start)', 'start')
row('Cap kept, 100 further iterations', 'kept')
row('Cap released, 100 further iterations', 'released')
cs = []
for c in conds:
    v = num[str(c)]
    cs.append(fmt_d(v['diff'], len(v['draws']), v['npos']) if 'diff' in v else '--')
rows.append('Released $-$ kept & ' + ' & '.join(cs) + '\\\\')
rows.append('\\midrule')
def fmt_pm(key, nkey):
    cs = []
    for c in conds:
        v = num[str(c)]
        cs.append('$%+.2f$ (%d/%d)' % (v[key][0], v[nkey], len(v['draws'])))
    return cs
rows.append('Released $-$ kept, objective $J$ (160 design samples) & ' + ' & '.join(fmt_pm('dJ', 'posJ')) + '\\\\')
rows.append('Released $-$ kept, rate on 2000 fresh samples & ' + ' & '.join(fmt_pm('dR', 'posR')) + '\\\\')
rows.append('$J-$fresh-sample rate: kept / released & ' + ' & '.join('$%+.2f$ / $%+.2f$' % (num[str(c)]['gap_kept'][0], num[str(c)]['gap_rel'][0]) for c in conds) + '\\\\')
rows.append('\\midrule')
row('Capped design from the closed form, 150 iterations', 'cold_cap')
row('Compression-aware from the closed form, 150 iterations', 'cold_aware')
open(S + 'new_tab_release_rows.tex', 'w').write('\n'.join(rows) + '\n')

# rho / pmax
R = collections.defaultdict(dict)
import glob
for f in glob.glob(S + 'rho_L*.log') + glob.glob(S + 'pmax*_L*.log'):
    for ln in open(f):
        m = re.match(r'RHO (\d+) (\d+) (-?\d+) ([\d.]+) (\d+) (\d+) ([\d.]+)', ln)
        if m: R[(int(m[1]), float(m[4]), int(m[5]))][int(m[2])] = float(m[7])
for L in (100, 200):
    for d in (1, 2, 3): R[(L, 0.1, 100)][d] = A.val('linreg', L, d, 0, 150)
ref_aware = {L: {d: A.val('aware', L, d, 0, 150) for d in (1, 2, 3)} for L in (100, 200)}
rr = []
def rrow(label, rho, pm):
    cs = []
    for L in (100, 200):
        dct = R.get((L, rho, pm), {})
        ds = [d for d in (1, 2, 3) if d in dct]
        cs.append(fmt(mean_se([dct[d] for d in ds])) + (' (n=%d)' % len(ds) if len(ds) != 3 else '') if ds else '--')
        num['rho_%g_pmax_%d_L%d' % (rho, pm, L)] = {'n': len(ds), 'mean_se': mean_se([dct[d] for d in ds]) if ds else None, 'per_draw': [dct[d] for d in ds]}
    rr.append(label + ' & ' + ' & '.join(cs) + '\\\\')
for rho in (0.03, 0.1, 0.3, 1.0): rrow('$\\rho_{\\max}=%g$, $p_{\\max}=100$' % rho, rho, 100)
rr.append('\\midrule')
for pm in (30, 100, 300): rrow('$\\rho_{\\max}=0.1$, $p_{\\max}=%d$' % pm, 0.1, pm)
cs = []
for L in (100, 200):
    ds = [d for d in (1, 2, 3) if d in ref_aware[L]]
    cs.append(fmt(mean_se([ref_aware[L][d] for d in ds])))
    num['aware_cold_L%d' % L] = mean_se([ref_aware[L][d] for d in ds])
rr.append('\\midrule')
rr.append('Compression-aware (no cap), $p_{\\max}=100$ & ' + ' & '.join(cs) + '\\\\')
open(S + 'new_tab_rho_rows.tex', 'w').write('\n'.join(rr) + '\n')
json.dump(num, open(S + 'new_numbers.json', 'w'), indent=1, default=str)
print(open(S + 'new_tab_release_rows.tex').read()); print(open(S + 'new_tab_rho_rows.tex').read())
