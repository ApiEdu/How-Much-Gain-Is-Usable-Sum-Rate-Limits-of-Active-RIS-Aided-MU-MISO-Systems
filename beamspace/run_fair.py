import numpy as np
from bsmodel import *
for K, L, Ns in [(4, 6, 6), (8, 16, 12)]:
    for fair in (0.0, 0.5):
        R, M, J = [], [], []
        for dr in range(4):
            sc = Scene(np.random.default_rng(100 + dr), K=K, L=L); sc.Pris = 0.1; A = np.sqrt(sc.Pris / sc.N)
            S = select_greedy(sc, Ns, np.inf, aware=True)
            d = design_refined(sc, S, A, fair=fair); s, _ = evaluate(sc, d['S'], d['W'], d['ph'], d['p'], A); r = np.log2(1 + s)
            R.append(r.sum()); M.append(r.min()); J.append(r.sum() ** 2 / (K * (r ** 2).sum()))
        print(f"K={K} L={L} {'sum-rate  ' if fair==0 else 'prop-fair '}: sum-rate {np.mean(R):5.2f}  min-user {np.mean(M):5.2f}  Jain {np.mean(J):.2f}", flush=True)
