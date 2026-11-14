# Beamspace (lens-array) MU system with an active RIS - practical study (work in progress)

Python (numpy/scipy). Model in `bsmodel.py`; experiments: `run_h1.py`, `diag_decouple.py`, `run_ee.py`, `run_ee2.py`, `run_k.py`.

## Status (14 Nov 2026) - preliminary, not for citation
- Beam-set ranking is the same with a linear and a saturating amplifier (Spearman >= 0.997): saturation-aware beam selection brings nothing in this model (`diag_decouple.py`).
- EE-optimal number of RF chains is about 1.5K for K <= 6 users; EE gain over full-digital 1.4x-3.9x depending on assumed power constants (`run_ee2.py`).
- Sum-rate does not scale beyond the effective rank of the BS-RIS channel G (L = 6 paths -> about 5 usable streams); with L = 16 K = 8 improves, K = 16 does not (`run_k.py`).

## Open points / assumptions
- Power constants are placeholders inside the literature ranges: RF chain 100-250 mW, active-RIS element DC 8-26.5 mW (see the paper draft for references); BS PA efficiency 0.4 and fixed power 1 W are assumptions.
- RIS phases by a heuristic alternating rule (no gradient/WMMSE refinement), uniform amplifier gain, RZF precoder, greedy beam selection, K <= 16, one geometry (28 GHz, 50 m / 20 m), 3-8 channel draws.
- Passive-RIS baseline not meaningful in this link budget (SNR about -6 dB); pick a geometry where it works before comparing.
