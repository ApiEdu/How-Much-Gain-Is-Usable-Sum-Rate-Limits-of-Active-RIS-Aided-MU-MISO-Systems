"""Diagnostic: do beam sets rank the same under a linear amplifier and under saturation?
Also: how uneven is the RIS input power Pu_n across elements (the only way saturation can interact with beam choice)."""
import numpy as np
from scipy.stats import spearmanr
import json
from bsmodel import *
OUT = []
for dr in range(4):
    sc = Scene(np.random.default_rng(100 + dr), dBR=20.0, dRU=10.0); rng = np.random.default_rng(dr)
    ph0 = phases(sc, np.arange(sc.M)); en = np.argsort(-cascaded_energy(sc, ph0))
    sets = [np.sort(rng.choice(en[:20], 6, replace=False)) for _ in range(120)]
    for beta in (-5, 0, 5):
        A = np.sqrt(10 ** (beta / 10) * sc.Pris / sc.N)
        lin, sat, spread = [], [], []
        for S in sets:
            d = design(sc, S, A, aware=True)
            sat.append(rate_analytic(sc, d, A))
            dl = design(sc, S, np.inf, aware=True)
            lin.append(rate_analytic(sc, dl, np.inf))
            Pu = (np.abs(sc.Gb[:, S] @ d['W']) ** 2).sum(1)
            spread.append(10 * np.log10(Pu.max() / Pu.min()))
        i_lin, i_sat = int(np.argmax(lin)), int(np.argmax(sat))
        OUT.append(dict(draw=dr, beta=beta, rho=float(spearmanr(lin, sat)[0]), same_best=bool(abs(sat[i_lin] - sat[i_sat]) < 1e-9), loss=float(sat[i_sat] - sat[i_lin]), spread=float(np.median(spread))))
        print(f"draw {dr} beta {beta:+d}: spearman(lin,sat)={spearmanr(lin,sat)[0]:.3f}  best-by-lin sat-rate={sat[i_lin]:.2f} vs best sat={sat[i_sat]:.2f}"
              f"  median Pu spread across elements={np.median(spread):.1f} dB", flush=True)

json.dump(OUT, open('decouple.json', 'w'))
print('min rho', min(o['rho'] for o in OUT), 'max loss', max(o['loss'] for o in OUT), 'same best', sum(o['same_best'] for o in OUT), '/', len(OUT))
