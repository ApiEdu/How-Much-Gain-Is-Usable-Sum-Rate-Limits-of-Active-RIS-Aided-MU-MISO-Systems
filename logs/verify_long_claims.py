"""Recompute every number quoted in the subsection "Longer Optimiser Budget" (tab:long) and in the sentences that cite it,
and assert it against the rounded value in the manuscript. Prints 'all long-budget claims verified' on success.
Run without -I:  python3 logs/verify_long_claims.py
"""
import os, re, glob
from decimal import Decimal, ROUND_HALF_UP
import analyze_long as A

HERE = os.path.dirname(os.path.abspath(__file__))


def r2(x):
    """Round half away from zero on the 3-decimal value (the logs have 3 decimals)."""
    return float(Decimal(f'{abs(x):.3f}').quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)) * (1 if x >= 0 else -1)


def close(a, b):
    assert abs(a - b) < 1e-9, (a, b)


data = A.load()
rows = A.main_quiet()
S = A.summary(rows)
pairs = [(L, d, b) for L in (100, 200) for b in (-5, 0) for d in (1, 2)]
assert len(pairs) == 8 and all((L, d, b, 'aware', 600) in data for (L, d, b) in pairs)
assert len(data) == 8 * 2 * 3                      # 8 pairs x 2 designs x 3 checkpoints

# --- 150-iteration checkpoints reproduce the earlier logs exactly (16 of 16 designs)
import re
ref = {}
for L in (100, 200):
    for fn, tag, des in ((f'convgrad_L{L}.log', 'CONV', 'aware'), (f'convlin_L{L}.log', 'CONVLIN', 'linreg')):
        for ln in open(os.path.join(HERE, fn)):
            m = re.match(tag + r' (\d+) (\d+) (-?\d+) it (\d+).*?simrate ([\d.]+)', ln)
            if m and int(m.group(4)) == 150:
                ref[(int(m.group(1)), int(m.group(2)), int(m.group(3)), des)] = float(m.group(5))
    for ln in open(os.path.join(HERE, f'main150_L{L}.log')):
        t = ln.split()
        if ln.startswith('M150 ') and t[5] == '150' and t[4] in ('aware', 'linreg'):
            ref[(int(t[1]), int(t[2]), int(t[3]), t[4])] = float(t[6])
n_match = 0
for (L, d, b) in pairs:
    for mode in ('aware', 'linreg'):
        assert (L, d, b, mode) in ref, (L, d, b, mode)
        assert abs(data[(L, d, b, mode, 150)][0] - ref[(L, d, b, mode)]) < 1e-9
        n_match += 1
assert n_match == 16

# --- sign of aware - cap at every checkpoint
delta = {(L, d, b, it): data[(L, d, b, 'aware', it)][0] - data[(L, d, b, 'linreg', it)][0]
         for (L, d, b) in pairs for it in A.CPS}
assert all(v > 0 for k, v in delta.items() if k[0] == 100)
assert all(v < 0 for k, v in delta.items() if k[0] == 200)

# --- two-draw means (tab:long and text)
for key, (a, c) in {(100, -5): (0.94, 0.82), (100, 0): (0.29, 0.20), (200, -5): (-0.16, -0.14), (200, 0): (-0.16, -0.13)}.items():
    assert r2(S[key + (150,)]['dr']) == a, (key, S[key + (150,)]['dr'])
    assert r2(S[key + (600,)]['dr']) == c, (key, S[key + (600,)]['dr'])
assert r2(S[(200, -5, 300)]['dr']) == -0.11 and r2(S[(200, 0, 300)]['dr']) == -0.14
assert r2(S[(100, -5, 300)]['dr']) == 0.92 and r2(S[(100, 0, 300)]['dr']) == 0.20

# --- table cells
tab = {(100, -5): (12.54, 12.76, 13.47, 13.58), (100, 0): (15.48, 15.78, 15.77, 15.98),
       (200, -5): (18.71, 19.01, 18.55, 18.87), (200, 0): (20.84, 21.23, 20.67, 21.10)}
for (L, b), (c150, c600, a150, a600) in tab.items():
    assert r2(S[(L, b, 150)]['cap']) == c150 and r2(S[(L, b, 600)]['cap']) == c600
    assert r2(S[(L, b, 150)]['aware']) == a150 and r2(S[(L, b, 600)]['aware']) == a600
per = {(100, -5): (0.78, 0.86), (100, 0): (0.16, 0.24), (200, -5): (-0.25, -0.03), (200, 0): (-0.21, -0.05)}
for (L, b), (p1, p2) in per.items():
    assert r2(delta[(L, 1, b, 600)]) == p1 and r2(delta[(L, 2, b, 600)]) == p2, (L, b)
dj = {(100, -5): 1.21, (100, 0): 0.65, (200, -5): 0.13, (200, 0): 0.04}
for k, v in dj.items():
    assert r2(S[k + (600,)]['dJ']) == v, (k, S[k + (600,)]['dJ'])

# --- change of the difference between 150 and 600 iterations, per pair
ch = [delta[(L, d, b, 600)] - delta[(L, d, b, 150)] for (L, d, b) in pairs]
assert r2(min(ch)) == -0.17 and r2(max(ch)) == 0.05, (min(ch), max(ch))
# shrinkage of the two-draw mean advantage at 100 m
assert r2(S[(100, -5, 150)]['dr'] - S[(100, -5, 600)]['dr']) == 0.12
assert r2(S[(100, 0, 150)]['dr'] - S[(100, 0, 600)]['dr']) == 0.09
# the cap gained more than the aware design at 100 m (mean over the two draws)
for b in (-5, 0):
    gc = S[(100, b, 600)]['cap'] - S[(100, b, 150)]['cap']
    ga = S[(100, b, 600)]['aware'] - S[(100, b, 150)]['aware']
    assert gc > ga
# range of the advantage at 600 iterations at 100 m
v100 = [delta[(100, d, b, 600)] for b in (-5, 0) for d in (1, 2)]
assert r2(min(v100)) == 0.16 and r2(max(v100)) == 0.86
# 200 m, 600 iterations: -0.25 to -0.03
v200 = [delta[(200, d, b, 600)] for b in (-5, 0) for d in (1, 2)]
assert r2(min(v200)) == -0.25 and r2(max(v200)) == -0.03

# --- both designs still improving
g = [data[(L, d, b, m, 600)][0] - data[(L, d, b, m, 150)][0] for (L, d, b) in pairs for m in ('aware', 'linreg')]
assert len(g) == 16 and r2(min(g)) == 0.07 and r2(max(g)) == 0.45, (min(g), max(g))
g2 = [data[(L, d, b, m, 600)][0] - data[(L, d, b, m, 300)][0] for (L, d, b) in pairs for m in ('aware', 'linreg')]
assert len(g2) == 16 and all(x > 0 for x in g2) and r2(min(g2)) == 0.02 and r2(max(g2)) == 0.17, (min(g2), max(g2))

# --- design-time objective at 600 iterations
dJ = {(L, d, b): data[(L, d, b, 'aware', 600)][1] - data[(L, d, b, 'linreg', 600)][1] for (L, d, b) in pairs}
assert sum(1 for v in dJ.values() if v > 0) == 7
neg = [k for k, v in dJ.items() if v <= 0]
assert neg == [(200, 1, 0)] and r2(dJ[(200, 1, 0)]) == -0.10
ex_a = [data[(L, d, b, 'aware', 600)][1] - data[(L, d, b, 'aware', 600)][0] for (L, d, b) in pairs]
ex_c = [data[(L, d, b, 'linreg', 600)][1] - data[(L, d, b, 'linreg', 600)][0] for (L, d, b) in pairs]
assert r2(min(ex_a)) == 0.15 and r2(max(ex_a)) == 0.69, (min(ex_a), max(ex_a))
assert r2(min(ex_c)) == -0.01 and r2(max(ex_c)) == 0.58, (min(ex_c), max(ex_c))
diff = [a - c for a, c in zip(ex_a, ex_c)]
assert all(x > 0 for x in diff) and r2(min(diff)) == 0.11 and r2(max(diff)) == 0.48, (min(diff), max(diff))
# J advantage larger than the score advantage in all eight pairs (same statement as 'diff > 0')
assert all(dJ[(L, d, b)] > delta[(L, d, b, 600)] for (L, d, b) in pairs)
# 200 m: aware higher J in 3 of 4 pairs, cap higher score in all 4
p200 = [(200, d, b) for b in (-5, 0) for d in (1, 2)]
assert sum(1 for k in p200 if dJ[k] > 0) == 3 and all(delta[k + (600,)] < 0 for k in p200)

print('all long-budget claims verified')
