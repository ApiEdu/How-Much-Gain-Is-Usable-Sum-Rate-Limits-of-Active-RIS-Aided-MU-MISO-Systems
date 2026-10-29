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
The analysis scripts contain an absolute path in the variable `S` (the build container); set it to your own directory before running them.

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

Wall-clock times quoted in the article were measured on a shared two-core container and are indicative only.

## Licence
MIT License (see the file `LICENSE`).
