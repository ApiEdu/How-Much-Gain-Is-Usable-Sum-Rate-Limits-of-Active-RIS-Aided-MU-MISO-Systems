"""User grouping / time-sharing when K exceeds the usable rank of the BS-RIS channel.
Users are split into groups of size g by a correlation-based greedy rule; groups are served in equal time slots, the RIS is
re-configured per slot.  Per-user average rate = R_k / (number of groups).  Beams per group: Ns = ceil(1.5 g) by linear greedy."""
import copy, sys, numpy as np
from bsmodel import *
K = int(sys.argv[1]); ND = int(sys.argv[2]); Ls = [int(x) for x in sys.argv[3].split(',')]
PW = dict(eta_bs=0.4, P_rf=0.25, P_fix=1.0, eta_ris=0.4, P_el_act=0.02)

def groups_by_correlation(f, g):
    K = f.shape[0]; fn = f / np.linalg.norm(f, axis=1, keepdims=True); Cc = np.abs(fn @ fn.conj().T)
    left = list(range(K)); out = []
    while left:
        grp = [left.pop(int(np.argmax(np.linalg.norm(f[left], axis=1))))]
        while len(grp) < g and left:
            j = int(np.argmin([Cc[i, grp].max() for i in left])); grp.append(left.pop(j))
        out.append(grp)
    return out

for L in Ls:
    print(f"\nK={K}, L={L}, {ND} draws, Pris=20 dBm   (time-averaged over the slots)")
    print("  group size g | slots | Ns/slot | avg sum-rate | avg min-user rate | avg power [W] | EE")
    for g in [4, 6, 8, K]:
        if g > K: continue
        res = []
        for dr in range(ND):
            sc = Scene(np.random.default_rng(100 + dr), K=K, L=L); sc.Pris = 0.1; A = np.sqrt(sc.Pris / sc.N)
            grps = groups_by_correlation(sc.f, g); G = len(grps)
            Rk = np.zeros(K); Ps = []
            for grp in grps:
                sub = copy.copy(sc); sub.K = len(grp); sub.f = sc.f[grp]
                Ns = int(np.ceil(1.5 * len(grp)))
                S = select_greedy(sub, Ns, np.inf, aware=True, cand=min(sc.M, max(24, 2 * Ns)))
                d = design_refined(sub, S, A); s, Pr = evaluate(sub, d['S'], d['W'], d['ph'], d['p'], A)
                Rk[grp] = np.log2(1 + s)
                Ps.append(sc.Pbs / PW['eta_bs'] + Ns * PW['P_rf'] + PW['P_fix'] + Pr / PW['eta_ris'] + sc.N * PW['P_el_act'])
            Rk = Rk / G
            res.append((Rk.sum(), Rk.min(), np.mean(Ps), Rk.sum() / np.mean(Ps), G, int(np.ceil(1.5 * min(g, K)))))
        a = np.array(res)
        print(f"  {g:12d} | {int(a[0,4]):5d} | {int(a[0,5]):7d} | {a[:,0].mean():12.2f} | {a[:,1].mean():17.3f} | {a[:,2].mean():13.2f} | {a[:,3].mean():.2f}", flush=True)
