"""Scaling with the number of users K: sum-rate and EE vs number of RF chains Ns (beamspace, active RIS) vs full-digital."""
import sys, copy, numpy as np
from bsmodel import *
Ks = [int(x) for x in sys.argv[1].split(",")]; ND = int(sys.argv[2]); pdbm = 20; LP = int(sys.argv[3]) if len(sys.argv) > 3 else 6
PW = dict(eta_bs=0.4, P_rf=0.25, P_fix=1.0, eta_ris=0.4, P_el_act=0.02)
for K in Ks:
    mult = [1.0, 1.5, 2.0, 3.0]
    NS = sorted({min(64, int(np.ceil(m * K))) for m in mult})
    out = {Ns: [] for Ns in NS + [64]}
    for dr in range(ND):
        sc = Scene(np.random.default_rng(100 + dr), K=K, L=LP); sc.Pris = 10 ** (pdbm / 10) * 1e-3
        A = np.sqrt(sc.Pris / sc.N)
        best = 0
        for Ns in NS:
            S = select_greedy(sc, Ns, np.inf, aware=True)
            d = design(sc, S, A, aware=True); s, Pr = evaluate(sc, S, d['W'], d['ph'], d['p'], A)
            R = np.log2(1 + s).sum(); best = max(best, R)
            P = sc.Pbs / PW['eta_bs'] + Ns * PW['P_rf'] + PW['P_fix'] + Pr / PW['eta_ris'] + sc.N * PW['P_el_act']
            out[Ns].append((R, P, float(np.min(np.log2(1 + s)))))
        allb = np.arange(64); d = design(sc, allb, A, aware=True); s, Pr = evaluate(sc, allb, d['W'], d['ph'], d['p'], A)
        R = max(best, np.log2(1 + s).sum())
        P = sc.Pbs / PW['eta_bs'] + 64 * PW['P_rf'] + PW['P_fix'] + Pr / PW['eta_ris'] + sc.N * PW['P_el_act']
        out[64].append((R, P, float(np.min(np.log2(1 + s)))))
    print(f"\nK={K} L={LP} (Pris={pdbm} dBm, {ND} draws)   Ns : sum-rate [b/s/Hz] | min-user rate | power [W] | EE")
    for Ns, v in out.items():
        R = np.mean([a for a, _, _ in v]); P = np.mean([b for _, b, _ in v]); mn = np.mean([c for _, _, c in v]); EE = np.mean([a / b for a, b, _ in v])
        print(f"  Ns={Ns:2d}: R={R:6.2f}  Rmin={mn:5.2f}  P={P:5.2f}  EE={EE:5.2f}", flush=True)
