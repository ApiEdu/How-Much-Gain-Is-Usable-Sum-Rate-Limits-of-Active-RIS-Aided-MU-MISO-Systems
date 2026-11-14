import json, sys, numpy as np, collections
f = sys.argv[1] if len(sys.argv) > 1 else 'results_S2.jsonl'
D = collections.defaultdict(list)
for ln in open(f):
    r = json.loads(ln); D[(r['K'], r['L'])].append(r['res'])
def ms(x): x = np.array(x); return x.mean(), x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0
for (K, L) in sorted(D):
    rs = D[(K, L)]; n = len(rs)
    print(f"\n=== K={K}, L={L}, n={n} draws (mean +- s.e.)  [scheme: sum-rate | min-user | Jain | power W | EE]")
    for k in rs[0]:
        v = lambda m: ms([r[k][m] for r in rs])
        s, mn, j, p, ee = v('sum'), v('mn'), v('jain'), v('P'), v('EE')
        print(f"  {k:14s} {s[0]:6.2f}+-{s[1]:4.2f} | {mn[0]:5.2f}+-{mn[1]:4.2f} | {j[0]:4.2f} | {p[0]:5.2f} | {ee[0]:5.2f}+-{ee[1]:4.2f}")
    act = [k for k in rs[0] if k.startswith('act')]
    best = max(act, key=lambda k: np.mean([r[k]['EE'] for r in rs]))
    full = [k for k in act if k.endswith('_n64')]; bf = max(full, key=lambda k: np.mean([r[k]['EE'] for r in rs]))
    ratio = [r[best]['EE'] / r[bf]['EE'] for r in rs]
    print(f"  -> best-EE active scheme: {best}; best full-digital-type scheme: {bf}; paired EE ratio {np.mean(ratio):.2f} (min {np.min(ratio):.2f}, max {np.max(ratio):.2f})")
