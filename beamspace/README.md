# Lens-antenna (beamspace) MU system aided by an active RIS: practical study

Python (numpy/scipy/matplotlib). Model and design routines: `bsmodel.py`. Manuscript draft: `paper/main.tex` (+ `paper/main.pdf`).

## Reproduce
```
python run_full.py results_S2.jsonl 100000 30 S2 4,8,16 6,16   # 180 tasks (30 draws x K in {4,8,16} x L in {6,16}), resumable, ~25 min on 2 cores
python validate_bussgang.py          # analytic (Bussgang) vs symbol-level rate  -> validation.json
python diag_decouple.py              # linear vs saturating beam-set ranking      -> decouple.json
python make_figs.py                  # figures, tables, paper/numbers.json
python build_paper.py                # fills paper/main.tex.in from the numbers and compiles paper/main.pdf
```
All numbers quoted in the manuscript are generated from `results_S2.jsonl`, `validation.json` and `decouple.json`; nothing is typed by hand.

## What the study is (and is not)
A simulation-based practical study, not a new algorithm: heuristic design (regularised ZF, gradient RIS-phase refinement on a proportional-fair objective, greedy beam selection on a linear-amplifier model, correlation-based user grouping with time sharing).
Limitations: K <= 16, one geometry (28 GHz, 20 m / 10 m, free space, direct link blocked), idealised lens (unitary DFT) and channels, perfect CSI, uniform amplifier gain, analytic (slightly optimistic) rates, several assumed power constants (eta_BS, eta_RIS, P_fix, P_el of a passive element), per-draw best N_s is an oracle choice.

## Earlier exploratory scripts (kept for traceability, superseded by `run_full.py`)
`run_h1.py`, `run_ee.py`, `run_ee2.py`, `run_k.py`, `run_refine.py`, `run_group.py`, `run_fair.py`, `analyze_*.py`.
The initial hypothesis (saturation-aware beam selection brings a large gain) was not supported: see `diag_decouple.py` and Section "Saturation and beam selection" of the manuscript.
