"""Option 1: energy efficiency of lens-based (beamspace) active-RIS systems -- RF chains vs RIS amplifier budget.
Rates: analytic Bussgang model (validated against symbol-level simulation to ~0.1 bit/s/Hz); beam sets chosen once per
(draw, Ns) by linear-model greedy selection (ranking is independent of saturation, see diag_decouple.py).
POWER MODEL = ASSUMPTIONS (only the 250 mW/RF-chain value comes from the literature; the rest are placeholders):
  P_total = Pbs/eta_bs + Ns*P_rf + P_fix + Pris_out/eta_ris + N*P_el
"""
import sys, copy, json, numpy as np
from bsmodel import *

ND = int(sys.argv[1]) if len(sys.argv) > 1 else 6
NS = [4, 6, 8, 12, 16]
PRIS_DBM = [-10, 0, 10, 20]
PW = dict(eta_bs=0.4, P_rf=0.25, P_fix=1.0, eta_ris=0.4, P_el_act=0.02, P_el_pas=0.001)

def total_power(Ns, Pout, active, pw, N, Pbs):
    P = Pbs / pw['eta_bs'] + Ns * pw['P_rf'] + pw['P_fix']
    P += (Pout / pw['eta_ris'] + N * pw['P_el_act']) if active else N * pw['P_el_pas']
    return P

res = {}   # (arch, Ns, Pris) -> list of (rate, power)
for dr in range(ND):
    base = Scene(np.random.default_rng(100 + dr))
    sets = {Ns: select_greedy(base, Ns, np.inf, aware=True) for Ns in NS}
    allb = np.arange(base.M)
    for pdbm in PRIS_DBM:
        sc = copy.copy(base); sc.Pris = 10 ** (pdbm / 10) * 1e-3
        A = np.sqrt(sc.Pris / sc.N)                      # beta = 0 dB
        for Ns in NS:
            d = design(sc, sets[Ns], A, aware=True)
            s, Pr = evaluate(sc, d['S'], d['W'], d['ph'], d['p'], A)
            res.setdefault(('active', Ns, pdbm), []).append((np.log2(1 + s).sum(), total_power(Ns, Pr, True, PW, sc.N, sc.Pbs)))
        # full-digital reference: M RF chains, may use any subset -> best rate over the tested subsets (and all beams)
        best = max(res[('active', Ns, pdbm)][-1][0] for Ns in NS)
        dA = design(sc, allb, A, aware=True); sA, PrA = evaluate(sc, allb, dA['W'], dA['ph'], dA['p'], A)
        rA = max(best, np.log2(1 + sA).sum())
        res.setdefault(('fulldig', sc.M, pdbm), []).append((rA, total_power(sc.M, PrA, True, PW, sc.N, sc.Pbs)))
    # passive RIS (no amplifier, no amplifier noise)
    sp = copy.copy(base); sp.sv2 = 0.0
    for Ns in NS:
        d = design(sp, sets[Ns], np.inf, aware=True, pgrid=np.array([1.0]))
        s, _ = evaluate(sp, d['S'], d['W'], d['ph'], 1.0, np.inf)
        res.setdefault(('passive', Ns, None), []).append((np.log2(1 + s).sum(), total_power(Ns, 0, False, PW, sp.N, sp.Pbs)))
    print('draw', dr, 'done', flush=True)

print('\nmean over', ND, 'draws.  R = sum-rate [bit/s/Hz], P = total power [W], EE = mean(R/P) [bit/s/Hz/W]')
rows = []
for (arch, Ns, pd), v in sorted(res.items(), key=lambda kv: (kv[0][0], str(kv[0][2]), kv[0][1])):
    R = np.mean([a for a, _ in v]); P = np.mean([b for _, b in v]); EE = np.mean([a / b for a, b in v])
    rows.append((arch, Ns, pd, R, P, EE))
    print(f"{arch:8s} Ns={Ns:2d} Pris={str(pd):>5s} dBm  R={R:6.2f}  P={P:5.2f}  EE={EE:5.2f}")
json.dump([list(map(lambda x: x if not isinstance(x, np.floating) else float(x), r)) for r in rows], open('ee_results.json', 'w'))
