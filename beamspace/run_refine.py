import copy, numpy as np
from bsmodel import *
for K, L, Ns in [(4, 6, 6), (8, 16, 12), (16, 16, 24)]:
    r0, r1, m0, m1 = [], [], [], []
    for dr in range(4):
        sc = Scene(np.random.default_rng(100 + dr), K=K, L=L); sc.Pris = 0.1; A = np.sqrt(sc.Pris / sc.N)
        S = select_greedy(sc, Ns, np.inf, aware=True)
        d0 = design(sc, S, A, aware=True); d1 = design_refined(sc, S, A)
        s0, _ = evaluate(sc, d0['S'], d0['W'], d0['ph'], d0['p'], A); s1, _ = evaluate(sc, d1['S'], d1['W'], d1['ph'], d1['p'], A)
        r0.append(np.log2(1 + s0).sum()); r1.append(np.log2(1 + s1).sum()); m0.append(np.log2(1 + s0).min()); m1.append(np.log2(1 + s1).min())
        sym = rate_symbol(sc, d1, A, T=20000, seed=3)
        print(f"K={K} L={L} draw {dr}: heuristic {r0[-1]:.2f} -> refined {r1[-1]:.2f} (symbol-level {sym:.2f}); min-user {m0[-1]:.2f} -> {m1[-1]:.2f}", flush=True)
    print(f"  MEAN K={K}: sum-rate {np.mean(r0):.2f} -> {np.mean(r1):.2f}; min-user {np.mean(m0):.2f} -> {np.mean(m1):.2f}")
