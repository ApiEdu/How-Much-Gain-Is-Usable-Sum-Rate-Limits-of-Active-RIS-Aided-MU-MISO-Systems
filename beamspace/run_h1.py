"""H1: does saturation-aware beam selection beat conventional energy-based selection under amplifier saturation?"""
import sys, numpy as np
from bsmodel import *

Ns_list = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [6]
betas = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [-5, 0, 5]
ndraw = int(sys.argv[3]) if len(sys.argv) > 3 else 3
for dr in range(ndraw):
    sc = Scene(np.random.default_rng(100 + dr))
    # sanity: linear-model SNR/rate with all beams
    full = np.arange(sc.M)
    for beta in betas:
        A = np.sqrt(10 ** (beta / 10) * sc.Pris / sc.N)
        dfull = design(sc, full, A, aware=True)
        row = dict(full=rate_symbol(sc, dfull, A, seed=1))
        for Ns in Ns_list:
            ph0 = phases(sc, full)
            Se = select_energy(sc, Ns, ph0)
            Sa = select_greedy(sc, Ns, A, aware=True)
            Sl = select_greedy(sc, Ns, A, aware=False)
            de = design(sc, Se, A, aware=True, ph=ph0)           # energy selection + aware gain
            da = design(sc, Sa, A, aware=True)                     # aware selection + aware gain
            dl = design(sc, Sl, A, aware=False)                    # linear-model selection + blind gain
            print(f"draw {dr} beta {beta:+.0f} Ns {Ns}  energy {rate_symbol(sc,de,A,seed=1):.2f}  aware {rate_symbol(sc,da,A,seed=1):.2f}"
                  f"  blind {rate_symbol(sc,dl,A,seed=1):.2f}  full(M={sc.M}) {row['full']:.2f}  "
                  f"| analytic-vs-symbol(aware) {rate_analytic(sc,da,A):.2f}/{rate_symbol(sc,da,A,seed=1):.2f}  p={da['p']:.1f}", flush=True)
