"""EE-ratio sensitivity to the assumed power constants (re-computed from the stored time-averaged power P)."""
import json, itertools, collections, numpy as np
D = collections.defaultdict(list)
for ln in open('results_S2.jsonl'):
    r = json.loads(ln); D[(r['K'], r['L'])].append(r['res'])
N = 128
print("EE gain of the best beamspace-active scheme over the full-digital (64 RF chains) scheme with the same grouping, mean over all (K,L,draw)")
print("  P_rf[mW]  P_el[mW] | gain   (min over draws)")
for Prf, Pel in itertools.product([0.1, 0.16, 0.25], [0.0085, 0.02, 0.0265]):
    gains = []
    for (K, L), rs in D.items():
        for r in rs:
            best = {}
            for k, v in r.items():
                if not k.startswith('act'): continue
                g = k.split('_')[1]; Ns = int(k.split('_n')[1])
                P = v['P'] + Ns * (Prf - 0.25) + N * (Pel - 0.02)
                best.setdefault(g, {})[Ns] = v['sum'] / P
            for g, d in best.items():
                if 64 in d and len(d) > 1: gains.append(max(v for n, v in d.items() if n != 64) / d[64])
    print(f"  {Prf*1e3:7.0f}  {Pel*1e3:8.1f} | x{np.mean(gains):.2f}  ({np.min(gains):.2f})")
