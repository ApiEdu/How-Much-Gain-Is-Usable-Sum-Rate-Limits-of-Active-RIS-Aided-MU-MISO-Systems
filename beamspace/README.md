# Beamspace (lens-array) MU system with an active RIS - practical study (work in progress)

Python (numpy/scipy). Model in `bsmodel.py`; experiments: `run_h1.py`, `diag_decouple.py`, `run_ee.py`, `run_ee2.py`, `run_k.py`.

## Status (14 Nov 2026) - preliminary, not for citation
- Beam-set ranking is the same with a linear and a saturating amplifier (Spearman >= 0.997): saturation-aware beam selection brings nothing in this model (`diag_decouple.py`).
- EE-optimal number of RF chains is about 1.5K for K <= 6 users; EE gain over full-digital 1.4x-3.9x depending on assumed power constants (`run_ee2.py`).
- Sum-rate does not scale beyond the effective rank of the BS-RIS channel G (L = 6 paths -> about 5 usable streams); with L = 16 K = 8 improves, K = 16 does not (`run_k.py`).

- Gradient refinement of the RIS phases (`refine_phases`, analytic gradient, W fixed) raises the sum-rate by 7-17% over the alternating heuristic; analytic and symbol-level rates agree within about 0.05 b/s/Hz (`run_refine.py`).
- Grouping users into slots of size g about equal to the usable rank (`run_group.py`, K = 16) raises the time-averaged sum-rate by 31-45% and EE by about 2x versus serving all K at once; the sum-rate objective still starves weak users (min-user rate about 0).
- A proportional-fair phase objective (`fair=0.5`) lifts the minimum user rate (0.02 -> 0.25 at K = 4; 0.61 -> 0.97 at K = 8) for a 1-4% sum-rate loss (`run_fair.py`); the precoder is still not fairness-aware.

## Open points / assumptions
- Power constants are placeholders inside the literature ranges: RF chain 100-250 mW, active-RIS element DC 8-26.5 mW (see the paper draft for references); BS PA efficiency 0.4 and fixed power 1 W are assumptions.
- RIS phases by a heuristic alternating rule (no gradient/WMMSE refinement), uniform amplifier gain, RZF precoder, greedy beam selection, K <= 16, one geometry (28 GHz, 50 m / 20 m), 3-8 channel draws.
- Passive-RIS baseline not meaningful in this link budget (SNR about -6 dB); pick a geometry where it works before comparing.
