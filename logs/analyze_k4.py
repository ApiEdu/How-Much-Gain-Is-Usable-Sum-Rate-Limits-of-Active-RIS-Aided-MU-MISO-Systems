"""K = 4 users, N = 256 (run_nk.m with NN=256 KK=4 GITERS=150).

Reads the Octave text-format result files res_k4a_L*_N256_K4.mat (draws 1-2) and
res_k4b_L*_N256_K4.mat (draws 3-4) from the repository root and prints, for each
distance L and saturation ratio beta, the mean +- standard error of the symbol-level
sum-rate of the four designs and the paired differences used in the paper.

Design index in R(draw, beta, design):
  1 aware   compression-aware gradient design, 150 iterations from the closed-form start
  2 linreg  linear-region-cap gradient design, 150 iterations from the capped closed form
  3 aware2  compression-aware design, 150 further iterations started from the linreg solution
            (300 iterations in total, so NOT an equal-budget comparison with 1 or 2)
  4 unaware linear-model gradient design (150 iterations), projected onto the feasible set
"""
import os, re, sys, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'   # repository root
NAMES = ['aware', 'linreg', 'aware2', 'unaware']
BETAS = [-5, 0, 5]


def read_octave_text_3d(path):
    """Minimal reader for 'type: matrix' arrays saved by Octave in text format (ndims 3)."""
    lines = open(path).read().split('\n')
    name = None
    for i, ln in enumerate(lines):
        if ln.startswith('# name: R'):
            name = i
            break
    assert name is not None, path
    j = name + 1
    while not lines[j].startswith('# ndims:'):
        j += 1
    dims = [int(x) for x in lines[j + 1].split()]
    vals = []
    j += 2
    while j < len(lines) and not lines[j].startswith('#') and lines[j].strip():
        vals.append(float(lines[j]))
        j += 1
    a = np.array(vals).reshape(dims, order='F')
    return a


def load(L):
    parts = []
    for tag in ('k4a', 'k4b'):
        f = ROOT + f'res_{tag}_L{L}_N256_K4.mat'
        if os.path.exists(f):
            parts.append(read_octave_text_3d(f))
    return np.concatenate(parts, axis=0) if parts else None   # (draws, beta, design)


def se(x):
    return x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else float('nan')


def main():
    out = {}
    for L in (100, 200):
        R = load(L)
        if R is None:
            continue
        n = R.shape[0]
        print(f'\n=== L = {L} m, K = 4, N = 256, n = {n} draws ===')
        for ib, b in enumerate(BETAS):
            x = R[:, ib, :]
            print(f' beta = {b:+d} dB')
            for s, nm in enumerate(NAMES):
                print(f'   {nm:8s} {x[:, s].mean():6.2f} +- {se(x[:, s]):.2f}   per-draw: ' +
                      ' '.join(f'{v:.2f}' for v in x[:, s]))
            pairs = [('linreg - unaware', 1, 3), ('aware - linreg (equal budget)', 0, 1),
                     ('aware2 - linreg (2x budget)', 2, 1), ('aware2 - aware (2x budget)', 2, 0)]
            for lab, a, c in pairs:
                d = x[:, a] - x[:, c]
                print(f'   {lab:32s} mean {d.mean():+.2f}  SE {se(d):.2f}  min {d.min():+.2f}  max {d.max():+.2f}  '
                      f'first higher in {int((d > 0).sum())}/{n}')
                out[(L, b, lab)] = d
    # overall ranges
    for lab in ('linreg - unaware', 'aware - linreg (equal budget)', 'aware2 - linreg (2x budget)'):
        v = np.concatenate([d for (L, b, l), d in out.items() if l == lab])
        print(f'\n{lab}: {len(v)} designs, range {v.min():+.2f} to {v.max():+.2f}, mean {v.mean():+.2f}, '
              f'positive in {int((v > 0).sum())}')


if __name__ == '__main__':
    main()
