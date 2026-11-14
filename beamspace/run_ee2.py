"""EE sensitivity: store raw (rate, RIS output power) once, then re-evaluate the (assumed) power model on a grid.
Rate model: analytic Bussgang with regularised ZF (alpha chosen per design); full-digital = best of subsets and all beams."""
import copy, itertools, numpy as np
from bsmodel import *
ND, NS, PD = 8, [4, 6, 8, 12, 16], [0, 10, 20, 25, 30, 35]
raw = {}
for dr in range(ND):
    base = Scene(np.random.default_rng(100 + dr))
    sets = {Ns: select_greedy(base, Ns, np.inf, aware=True) for Ns in NS}
    for pdbm in PD:
        sc = copy.copy(base); sc.Pris = 10 ** (pdbm / 10) * 1e-3; A = np.sqrt(sc.Pris / sc.N)
        rs = {}
        for Ns in NS:
            d = design(sc, sets[Ns], A, aware=True); s, Pr = evaluate(sc, d['S'], d['W'], d['ph'], d['p'], A)
            rs[Ns] = (np.log2(1 + s).sum(), Pr)
        allb = np.arange(sc.M); d = design(sc, allb, A, aware=True); s, Pr = evaluate(sc, allb, d['W'], d['ph'], d['p'], A)
        rs[64] = (max(np.log2(1 + s).sum(), *(v[0] for v in rs.values())), Pr)
        for Ns, v in rs.items(): raw.setdefault((Ns, pdbm), []).append(v)
    print('draw', dr, flush=True)
N, Pbs = 128, 1.0
def EE(Ns, pdbm, Pel, Prf, eta_bs=0.4, eta_ris=0.4, Pfix=1.0):
    v = raw[(Ns, pdbm)]
    return np.mean([R / (Pbs / eta_bs + Ns * Prf + Pfix + Pout / eta_ris + N * Pel) for R, Pout in v])
print("\nmean sum-rate [bit/s/Hz] (rows Ns, cols Pris dBm):", PD)
for Ns in NS + [64]: print(f"Ns={Ns:2d}", ' '.join(f"{np.mean([a for a,_ in raw[(Ns,p)]]):6.2f}" for p in PD))
print("\nsensitivity of the EE-optimal design (Ns*, Pris*) and its EE gain over full-digital (best Pris)")
print(" P_el[mW] P_rf[mW] |  Ns*  Pris*[dBm]  EE*   EEfull  gain")
for Pel, Prf in itertools.product([0.005, 0.02, 0.05], [0.1, 0.25, 0.4]):
    best = max(((EE(Ns, p, Pel, Prf), Ns, p) for Ns in NS for p in PD))
    full = max(EE(64, p, Pel, Prf) for p in PD)
    print(f" {Pel*1e3:7.0f} {Prf*1e3:8.0f} | {best[1]:4d} {best[2]:8d}   {best[0]:5.2f}  {full:5.2f}  x{best[0]/full:.2f}")
