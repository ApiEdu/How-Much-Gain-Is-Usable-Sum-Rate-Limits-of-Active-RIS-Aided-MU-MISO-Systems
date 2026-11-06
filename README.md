# How Much Gain Is Usable? Sum-Rate Limits of Active RIS-Aided MU-MISO Systems Under Amplifier Saturation

Code, simulation logs and analysis scripts for the article by A. B. Ibrahim and I. K. Eltahir (submitted to *Electronics*, MDPI).

## Contents
- Repository root: GNU Octave / MATLAB-compatible code (channel model, Rapp + Bussgang compression model, closed-form and gradient designs, genetic algorithm, symbol-level scoring; `*.m`) and Python code (WMMSE reference `fp_j98.py`, WMMSE with the gain cap inside the loop `fp_cap.py`, analysis scripts `analyze*.py`), plus the saved result files `res_*.mat`.
- `ch/` : saved channel realisations and WMMSE designs (`ch_L<L>_d<draw>.mat`, `fp*_L<L>_d<draw>.mat`) read by the evaluation scripts.
- `v6/` : the `res_*.mat` result files converted to MATLAB v6 format so that SciPy can read them (used by `analyze.py`, `analyze6k.py`).
- `logs/` : raw run logs and the scripts that turn them into the numbers, tables and figures of the article.

The scripts use relative paths (`ch/`, `v6/`): run them from the repository root.

## Requirements
GNU Octave 8.4.0 (or MATLAB); Python 3 with numpy, scipy, matplotlib and cvxpy (Clarabel solver).
The older analysis scripts (`analyze_conv.py`, `analyze_nk.py`, `analyze_val.py`, `make_*.py`) still contain an absolute path in the variable `S` (the build container); set it to your own directory before running them.
The scripts added on 4 November 2026 (listed below) and `logs/analyze_main150.py` find the logs from their own location, so they run from any directory. In `analyze_main150.py` this was the only change: the paths of the logs and of `v6/res_fp6k_L*.mat` are now derived from the file location; the numbers it prints are unchanged.

## Which log produces which result (table/figure numbers of the MDPI version)
| Result | Logs / scripts |
|---|---|
| Table 5, model validation | `valid_L*.log`, `analyze_val.py` |
| Tables 6, 7 and Figure 1 (14 channel draws) | draws 1-4: `study_L*.log`, `extra_L*.log`; draws 5-14: `more_L*.log`, `more2_L*.log`, `more3_L*.log`; WMMSE: `fp6k_L*.log`; combined by `analyze6k.py` (output `logs/summary_final6k.txt`); `make_fig_main_mdpi.py` redraws Figure 1 from the numbers in Tables 6 and 7 of the paper source |
| Tables 8, 9 (150 iterations) | `main150_L*.log`, `analyze_main150.py`, `make_main150_tex.py`, `analyze_conv.py` |
| In-loop WMMSE rows of Table 8 | `fpcap_*.log`, `fpcap_scores.txt` |
| Table 10 (N = 128) | `convnk_*.log`, `analyze_nk.py` |
| Tables 11, 12 (channel-estimation error) | `csi_L*.log`, `csig_L*.log`, `make_csi_tex.py` |
| Table 13 (genetic algorithm) | `study_L*.log` |
| Figure 2 (beta up to +30 dB) | `hbcf_L*.log`, `hbg_L*.log`, `analyze_hibeta.py`, `make_fig_hibeta_mdpi.py` |
| Grid check at beta = +30 dB | `gridtest.log` |

## Added on 4 November 2026: release-the-cap test, parameter sweeps, equal-budget draws 5-10, user separation
These files support the sections added to the revised manuscript (section titles are given because the table numbers differ between the journal templates).

New Octave code (repository root; run from the repository root, `octave-cli --no-gui <file>.m`, settings are read from environment variables):
- `grad_design_w.m` : the same Adam / finite-difference optimiser as `grad_design_h.m`, but the starting design is itself a candidate, so a run can be continued from an earlier design (used by the release-the-cap test).
- `run_warm.m` : release-the-cap test. A capped design is optimised for `CITERS` iterations; two continuations of `XITERS` iterations follow from that design, one with the gain cap kept and one with it released (compression-aware). Environment: `LV` (100 or 200), `CITERS=50`, `XITERS=100`, `DRAWS` (for example `1:4`), `BETAS` (for example `"-5 0 5"`). It prints `WARM` (symbol-level sum-rate), `PST` (gain statistics) and `DJ` (design-time objective).
- `run_rho.m` : sensitivity of the capped design to `rho_max` (environment `RHOS`) and to `p_max` (environment `PMAX`); `LV`, `GITERS=150`, `DRAWS`, `BETAS`.
- `run_main150x.m` : the equal-budget run of `run_main150.m` for further draws (`DRAWS=5:7`, `GITERS=150`, `LV`).
- `run_sep.m` : user-separation test (second user at (L, -SEP, 0); `SEP=30`, `GITERS=50`, `DRAWS=1:6`, `LV`). Its header comment is inherited from `run_main150x.m`.
- `run_all.sh L` (the commands are those that produced the logs; only the paths were changed to the repository layout): runs the release test and then the `rho_max` / `p_max` sweeps for one distance. The first run was stopped by a time limit of the compute environment after most of the sweeps; `run_rest.sh L` ran what was missing (`rho_max` = 1 for draw 3, `p_max` = 30 and 300).
- `GITERS`, `CITERS` and `XITERS` should be multiples of 10 (and `XITERS` a multiple of 20): the optimiser stores the design every 10 iterations and the scripts score the designs at `XITERS/2` and `XITERS`.

New analysis scripts (`logs/`): `analyze_warm.py`, `analyze_dj.py`, `stats_release.py`, `stats_rho.py` (release test and sweeps), `gen_tables_new.py` (writes the table rows `new_tab_release_rows.tex`, `new_tab_rho_rows.tex` and `new_numbers.json`), `verify_new_claims.py` (recomputes every number quoted in the new sections from the logs and checks it against the rounded value in the text; prints `ALL CLAIMS VERIFIED`), `stats10.py` and `analyze14.py` (equal-budget statistics over draws 1-10; `analyze14.py` is written for up to 14 draws, but 150-iteration logs exist only for draws 1-10, so it reports n = 10). Run them from anywhere, for example `python3 logs/verify_new_claims.py`.

| Result (section of the revised manuscript) | Logs / scripts |
|---|---|
| Equal-budget comparison, draws 1 to 10 (150 iterations) | draws 1-4: `main150_L*.log`; draws 5-10: `ext150_L*_d5-7.log`, `ext150_L*_d8-10.log`; `stats10.py`, `analyze14.py` |
| Releasing the cap from a common start | `warm_L*.log` (24 designs per distance: 4 draws x 3 beta x 2 continuations), `analyze_warm.py`, `analyze_dj.py`, `stats_release.py`, `gen_tables_new.py` |
| Sensitivity to the cap parameters (`rho_max`, `p_max`) | `rho_L*.log`, `rho_L*_b.log`, `pmax30_L*.log`, `pmax300_L*.log` (beta = 0 dB, draws 1-3; the reference values `rho_max` = 0.1, `p_max` = 100 come from `main150_L*.log`), `stats_rho.py`, `gen_tables_new.py` |
| Sensitivity to user separation (30 m) | `sep30_L*.log` (draws 1-6) |
| Check of the quoted numbers | `verify_new_claims.py` |

The logs begin with Octave warnings that `conv.m` and `conv2.m` of this repository shadow built-in functions; they do not affect the results.

Wall-clock times quoted in the article were measured on a shared two-core container and are indicative only.

## Added on 6 November 2026: four users (K = 4) at N = 256
These files support the subsection "Four Users at N = 256" of the revised manuscript.

- `run_nk.m` (already in the repository, unchanged) was run with `NN=256 KK=4 GITERS=150 BDB="-5 0 5"`, `LSET=100` or `200`, and `D0`/`D1` selecting the draws (`D0=0 D1=2` gives draws 1-2, `D0=2 D1=4` draws 3-4), `TAG=k4a` or `k4b`. The four users are at (L, 0, 0), (L+1, 0, 0), (L+2, 0, 0) and (L+3, 0, 0) m. Designs in the output: `aware` (compression-aware, 150 iterations from the closed-form start), `linreg` (capped, 150 iterations from the capped closed form), `aware2` (compression-aware, 150 further iterations started from the capped solution, 300 in total, so not an equal-budget comparison) and `unaware` (linear model, 150 iterations, projected onto the feasible set). Each draw takes about 40 minutes per distance on one core of a shared container (indicative).
- Result files (Octave text format, repository root): `res_k4a_L100_N256_K4.mat`, `res_k4a_L200_N256_K4.mat`, `res_k4b_L100_N256_K4.mat`, `res_k4b_L200_N256_K4.mat`. Logs: `logs/k4_L100_d1-2.log`, `logs/k4_L200_d1-2.log`, `logs/k4_L100_d3-4.log`, `logs/k4_L200_d3-4.log` (the Octave warnings at the top are the harmless `conv.m` / `conv2.m` shadowing warnings).
- Scripts (`logs/`): `analyze_k4.py` (mean +- standard error and paired differences per condition; reads the `.mat` files), `make_k4_tex.py` (prints the rows of the manuscript table), `verify_k4_claims.py` (recomputes every number quoted in the subsection and asserts it against the rounded value in the text; prints `all K=4 claims verified`). Run, for example, `python3 logs/verify_k4_claims.py`.

## Added on 6 November 2026: long-budget check (600 iterations)
These files support the subsection on the long-budget check of the revised manuscript (K = 2, N = 256, no other change to the model).

- `run_long.m` (repository root) runs the capped (`linreg`) and the compression-aware (`aware`) gradient design for 600 iterations, each from its own closed-form start as in `run_main150.m`, for one channel draw `DRAW` at distance `LV`, for the saturation ratios in `BETAS` (dB), and scores the iterate at the checkpoints in `CPS` (multiples of 10) by symbol-level simulation. Example: `LV=100 DRAW=1 GITERS=600 BETAS="-5 0" CPS="150 300 600" octave-cli --no-gui run_long.m > logs/long_L100_d1.log`. Output lines: `LONG <L> <draw> <beta_dB> <mode> <iteration> <symbol-level sum-rate> <design-time objective J>`. One run takes roughly an hour on one core of a shared container (indicative).
- Logs: `logs/long_L100_d1.log`, `logs/long_L200_d1.log`, `logs/long_L100_d2.log`, `logs/long_L200_d2.log` (`d1`, `d2` = channel draw).
- Scripts (`logs/`): `analyze_long.py` prints, for every distance, draw and saturation ratio, the rates after 150, 300 and 600 iterations, the aware-minus-cap differences and the design-time objective J; `make_long_tex.py` prints the rows of the manuscript table; `verify_long_claims.py` recomputes every number quoted in the subsection (including that the 150-iteration checkpoints reproduce the earlier 150-iteration logs exactly) and prints `all long-budget claims verified`. Run, for example, `cd logs; python3 verify_long_claims.py` (without `-I`, because the scripts import each other).

## Licence
MIT License (see the file `LICENSE`).
