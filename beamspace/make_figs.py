"""Figures, tables and the numbers quoted in the paper -> paper/figs/*.pdf, paper/numbers.json, paper/tab_*.tex"""
import json, collections, itertools, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from bsmodel import Scene
os.makedirs('paper/figs', exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.grid': True, 'grid.alpha': .3, 'figure.dpi': 150, 'lines.markersize': 4})
D = collections.defaultdict(list)
for ln in open('results_S2.jsonl'):
    r = json.loads(ln); D[(r['K'], r['L'])].append(r['res'])
n = min(len(v) for v in D.values()); 
for k in D: D[k] = D[k][:n]
def ms(x): x = np.array(x, float); return x.mean(), (x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0)
def stat(K, L, key, m): return ms([r[key][m] for r in D[(K, L)]])
NUM = dict(n=n)
# ---- Fig 1: singular values of G
fig, ax = plt.subplots(figsize=(3.4, 2.3))
for L, c in ((6, 'C0'), (16, 'C1')):
    sv = []
    for dr in range(20):
        sc = Scene(np.random.default_rng(100 + dr), L=L, dBR=20.0, dRU=10.0); s = np.linalg.svd(sc.G, compute_uv=False); sv.append(20 * np.log10(s[:24] / s[0] + 1e-18))
    sv = np.array(sv); ax.plot(np.arange(1, 25), np.median(sv, 0), 'o-', color=c, label=f'$L={L}$ paths (median of 20)')
    NUM[f'effrank20_L{L}'] = float(np.median([(row > -20).sum() for row in sv]))
ax.set_ylim(-60, 2); ax.set_xlabel('singular value index'); ax.set_ylabel('$20\\log_{10}(\\sigma_i/\\sigma_1)$ [dB]'); ax.legend(); fig.tight_layout(); fig.savefig('paper/figs/fig_svd.pdf'); plt.close(fig)
# ---- Fig 2: group size effect K=16 (Ns = 1.5 g) : sum-rate and min-user rate
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.3))
for L, c in ((6, 'C0'), (16, 'C1')):
    gs = [4, 8, 16]; keys = ['act_g4_n6', 'act_g8_n12', 'act_g16_n24']
    for ax, m, yl in ((axs[0], 'sum', 'avg. sum-rate [bit/s/Hz]'), (axs[1], 'mn', 'avg. min-user rate [bit/s/Hz]')):
        v = [stat(16, L, k, m) for k in keys]
        ax.errorbar(gs, [a for a, _ in v], yerr=[1.96 * b for _, b in v], fmt='o-', color=c, capsize=2, label=f'$L={L}$'); ax.set_xlabel('group size $g$ ($K=16$, $N_s=1.5g$)'); ax.set_ylabel(yl); ax.set_xticks(gs)
    for g, k in zip(gs, keys): NUM[f'K16_L{L}_g{g}_sum'] = stat(16, L, k, 'sum'); NUM[f'K16_L{L}_g{g}_mn'] = stat(16, L, k, 'mn'); NUM[f'K16_L{L}_g{g}_jain'] = stat(16, L, k, 'jain')
axs[0].legend(); fig.tight_layout(); fig.savefig('paper/figs/fig_group.pdf'); plt.close(fig)
# ---- Fig 3: EE vs number of RF chains (g = 4)
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.3)); Ns_l = [6, 12, 24, 64]
for (K, L), mk in zip([(8, 6), (8, 16), (16, 6), (16, 16)], ['o-', 's-', '^--', 'v--']):
    v = [stat(K, L, f'act_g4_n{N}', 'EE') for N in Ns_l]; axs[0].errorbar(range(4), [a for a, _ in v], yerr=[1.96 * b for _, b in v], fmt=mk, capsize=2, label=f'$K={K},L={L}$')
    v = [stat(K, L, f'act_g4_n{N}', 'sum') for N in Ns_l]; axs[1].errorbar(range(4), [a for a, _ in v], yerr=[1.96 * b for _, b in v], fmt=mk, capsize=2, label=f'$K={K},L={L}$')
for ax, yl in ((axs[0], 'EE [bit/s/Hz/W]'), (axs[1], 'avg. sum-rate [bit/s/Hz]')): ax.set_xticks(range(4)); ax.set_xticklabels(Ns_l); ax.set_xlabel('RF chains $N_s$ (64 = full digital), $g=4$'); ax.set_ylabel(yl)
axs[0].legend(fontsize=6); fig.tight_layout(); fig.savefig('paper/figs/fig_ee.pdf'); plt.close(fig)
# ---- paired EE ratio best-beamspace / full-digital (same grouping), per (K,L) and overall; rate loss
ratios, losses, per = [], [], {}
for (K, L), rs in sorted(D.items()):
    rr = []
    for r in rs:
        by = collections.defaultdict(dict)
        for k, v in r.items():
            if k.startswith('act'): by[k.split('_')[1]][int(k.split('_n')[1])] = v
        for g, d in by.items():
            if 64 in d and len(d) > 1:
                bn = max((nn for nn in d if nn != 64), key=lambda nn: d[nn]['EE']); rr.append(d[bn]['EE'] / d[64]['EE']); losses.append(1 - d[bn]['sum'] / d[64]['sum'])
    per[(K, L)] = ms(rr); ratios += rr
NUM['ee_ratio_all'] = ms(ratios); NUM['ee_ratio_min'] = float(np.min(ratios)); NUM['ee_ratio_max'] = float(np.max(ratios)); NUM['rate_loss_mean'] = float(np.mean(losses)); NUM['rate_loss_max'] = float(np.max(losses)); NUM['n_pairs'] = len(ratios)
# ---- tables
def f(x, d=2): return f'{x:.{d}f}'
with open('paper/tab_main.tex', 'w') as fo:
    for L in (6, 16):
        for k, lab in (('act_g4_n6', '$g{=}4,N_s{=}6$'), ('act_g4_n12', '$g{=}4,N_s{=}12$'), ('act_g4_n64', '$g{=}4,N_s{=}64$'), ('act_g8_n12', '$g{=}8,N_s{=}12$'), ('act_g8_n64', '$g{=}8,N_s{=}64$'), ('act_g16_n24', '$g{=}16,N_s{=}24$'), ('act_g16_n64', '$g{=}16,N_s{=}64$'), ('pas_g4_n12', 'passive RIS, $g{=}4,N_s{=}12$')):
            S, Mn, J, P, E = (stat(16, L, k, m) for m in ('sum', 'mn', 'jain', 'P', 'EE'))
            fo.write(f"{L} & {lab} & {f(S[0],1)}$\\pm${f(1.96*S[1],1)} & {f(Mn[0])} & {f(J[0])} & {f(P[0],1)} & {f(E[0])}$\\pm${f(1.96*E[1],2)} \\\\\n")
        if L == 6: fo.write('\\midrule\n')
with open('paper/tab_ratio.tex', 'w') as fo:
    for (K, L), (m, s) in sorted(per.items()): fo.write(f"{K} & {L} & {f(m)}$\\pm${f(1.96*s)} \\\\\n")
# sensitivity
N = 128; sens = []
for Prf, Pel in itertools.product([0.1, 0.16, 0.25], [0.0085, 0.02, 0.0265]):
    g = []
    for rs in D.values():
        for r in rs:
            best = {}
            for k, v in r.items():
                if k.startswith('act'): best.setdefault(k.split('_')[1], {})[int(k.split('_n')[1])] = v['sum'] / (v['P'] + int(k.split('_n')[1]) * (Prf - 0.25) + N * (Pel - 0.02))
            for gg, d in best.items():
                if 64 in d and len(d) > 1: g.append(max(v for nn, v in d.items() if nn != 64) / d[64])
    sens.append((Prf, Pel, float(np.mean(g)), float(np.min(g))))
with open('paper/tab_sens.tex', 'w') as fo:
    for a, b, m, mn in sens: fo.write(f"{a*1e3:.0f} & {b*1e3:.1f} & {m:.2f} & {mn:.2f} \\\\\n")
for _t in ('tab_main', 'tab_ratio', 'tab_sens'):   # last row without terminator (the template adds it): avoids "Misplaced \noalign"
    _x = open(f'paper/{_t}.tex').read().rstrip()
    assert _x.endswith('\\\\'); open(f'paper/{_t}.tex', 'w').write(_x[:-2].rstrip() + '\n')
NUM['sens_min_mean'] = min(s[2] for s in sens); NUM['sens_max_mean'] = max(s[2] for s in sens); NUM['sens_min_min'] = min(s[3] for s in sens)
# passive vs active (K=8,L=16, g=4, Ns=12)
NUM['pas_sum_K8L16'] = stat(8, 16, 'pas_g4_n12', 'sum'); NUM['act_sum_K8L16'] = stat(8, 16, 'act_g4_n12', 'sum')
NUM['pas_sum_max'] = max(max(r['pas_g4_n12']['sum'] for r in rs) for rs in D.values())
NUM['act_sum_min12'] = min(min(r['act_g4_n12']['sum'] for r in rs) for rs in D.values())
NUM['best_g_K16'] = {L: max([4, 8, 16], key=lambda g: stat(16, L, {4: 'act_g4_n6', 8: 'act_g8_n12', 16: 'act_g16_n24'}[g], 'sum')[0]) for L in (6, 16)}
NUM['K8_L16'] = {k: stat(8, 16, k, 'sum') for k in ('act_g4_n12', 'act_g8_n12', 'act_g8_n24')}
NUM['K8_L6'] = {k: stat(8, 6, k, 'sum') for k in ('act_g4_n12', 'act_g8_n12')}
NUM['ee_best'] = {f'{K}_{L}': max(((stat(K, L, k, 'EE')[0], k) for k in D[(K, L)][0] if k.startswith('act')))  for (K, L) in D}
# best configuration per group size at K = 16 (by sum-rate and by EE)
for L in (6, 16):
    for g in (4, 8, 16):
        ks = [k for k in D[(16, L)][0] if k.startswith(f'act_g{g}_') and not k.endswith('_n64')] + [f'act_g{g}_n64']
        bs = max(ks, key=lambda k: stat(16, L, k, 'sum')[0]); be = max(ks, key=lambda k: stat(16, L, k, 'EE')[0])
        NUM[f'K16_L{L}_g{g}_bestsum'] = [*stat(16, L, bs, 'sum'), bs]; NUM[f'K16_L{L}_g{g}_bestee'] = [*stat(16, L, be, 'EE'), be]
# all K=16 configurations (for text quoting)
NUM['K16_all'] = {L: {k: {m: stat(16, L, k, m) for m in ('sum', 'mn', 'jain', 'P', 'EE')} for k in D[(16, L)][0]} for L in (6, 16)}
# fixed-configuration paired EE ratios (no per-draw oracle choice of N_s)
fx = {}
for g, Ns_ in [(4, 6), (4, 12), (4, 24), (8, 12), (8, 24)]:
    rr = []
    for (K, L), rs in D.items():
        for r in rs:
            a, b = r.get(f'act_g{g}_n{Ns_}'), r.get(f'act_g{g}_n64')
            if a and b: rr.append(a['EE'] / b['EE'])
    fx[f'g{g}_n{Ns_}'] = (ms(rr), float(np.min(rr)), len(rr))
NUM['fixed_ratio'] = fx
# runtime per scene (all configurations of one (K,L,draw) task), shared 2-core CPU
sec = collections.defaultdict(list)
for ln in open('results_S2.jsonl'):
    r = json.loads(ln); sec[r['K']].append(r['sec'])
NUM['sec_by_K'] = {K: [float(np.median(v)), float(np.max(v))] for K, v in sorted(sec.items())}
json.dump(NUM, open('paper/numbers.json', 'w'), indent=1, default=float)
print(json.dumps(NUM, indent=1, default=float)[:3000])
