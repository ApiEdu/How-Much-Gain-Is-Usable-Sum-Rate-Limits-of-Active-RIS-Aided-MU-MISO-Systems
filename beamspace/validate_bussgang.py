"""Analytic (Bussgang) sum-rate versus symbol-level simulation (Gaussian symbols through the true Rapp amplifier)."""
import json, numpy as np
from bsmodel import *
rows = []
for K, L, Ns in [(4, 6, 12), (8, 16, 12), (16, 16, 24)]:
    for dr in range(10):
        sc = Scene(np.random.default_rng(100 + dr), K=K, L=L, dBR=20.0, dRU=10.0); sc.Pris = 0.1; A = np.sqrt(sc.Pris / sc.N)
        S = select_greedy(sc, Ns, np.inf, aware=True)
        d = design_refined(sc, S, A, fair=0.5, outer=2, steps=40)
        ra = rate_analytic(sc, d, A); rs = rate_symbol(sc, d, A, T=60000, seed=7 + dr)
        rows.append((K, L, ra, rs))
a = np.array(rows)
diff = a[:, 2] - a[:, 3]
out = dict(n=len(a), mean_abs=float(np.abs(diff).mean()), max_abs=float(np.abs(diff).max()), mean_rel=float(np.mean(np.abs(diff) / a[:, 3])), max_rel=float(np.max(np.abs(diff) / a[:, 3])), mean_signed=float(diff.mean()))
print(out); json.dump(out, open('validation.json', 'w'))
