#!/usr/bin/env python3
"""Recompute from the K = 4 result files every number quoted in the manuscript's 'Four Users' subsection and
assert that it matches the rounded value written there. Prints OK per claim; raises on mismatch."""
import analyze_k4 as A
import numpy as np

R = {L: A.load(L) for L in (100, 200)}
assert all(R[L].shape == (4, 3, 4) for L in R)          # 4 draws x 3 beta x 4 designs
AW, LR, A2, UN = 0, 1, 2, 3
def chk(name, val, expect, tol=0.005):
    ok = abs(val - expect) <= tol
    print(('OK  ' if ok else 'FAIL'), name, round(float(val), 4), 'expected', expect)
    assert ok, name
def d(L, ib, a, c): return R[L][:, ib, a] - R[L][:, ib, c]
allpairs = lambda a, c: np.concatenate([d(L, ib, a, c) for L in (100, 200) for ib in range(3)])

cl = allpairs(LR, UN)
chk('cap higher than linear model, all 24', (cl > 0).sum(), 24, 0)
chk('cap-linear min at 100 m', min(np.concatenate([d(100, ib, LR, UN) for ib in range(3)])), 1.37)
chk('cap-linear max at 100 m', max(np.concatenate([d(100, ib, LR, UN) for ib in range(3)])), 4.61)
chk('cap-linear min at 200 m', min(np.concatenate([d(200, ib, LR, UN) for ib in range(3)])), 13.88)
chk('cap-linear max at 200 m', max(np.concatenate([d(200, ib, LR, UN) for ib in range(3)])), 20.00)
for L, ib, au, cu, sh in [(100, 0, 5.74, 1.65, 29), (100, 1, 6.91, 4.09, 59), (100, 2, 5.38, 4.24, 79),
                          (200, 0, 16.86, 15.99, 95), (200, 1, 17.41, 17.75, 102), (200, 2, 14.26, 15.01, 105)]:
    a_, c_ = d(L, ib, AW, UN).mean(), d(L, ib, LR, UN).mean()
    chk(f'aware-linear L={L} b#{ib}', a_, au); chk(f'cap-linear L={L} b#{ib}', c_, cu)
    chk(f'share of gap L={L} b#{ib} (%)', 100 * c_ / a_, sh, 0.5)
ac = allpairs(AW, LR)
chk('aware-cap min', ac.min(), -1.22); chk('aware-cap max', ac.max(), 4.47); chk('aware-cap mean', ac.mean(), 1.30)
chk('aware-cap positive', (ac > 0).sum(), 17, 0)
a100 = np.concatenate([d(100, ib, AW, LR) for ib in range(3)])
chk('aware-cap at 100 m min', a100.min(), 0.98); chk('aware-cap at 100 m max', a100.max(), 4.47); chk('aware-cap at 100 m positive', (a100 > 0).sum(), 12, 0)
for L, ib, mean, wins in [(100, 0, 4.09, 4), (100, 1, 2.82, 4), (100, 2, 1.15, 4), (200, 0, 0.86, 4), (200, 1, -0.34, 1), (200, 2, -0.76, 0)]:
    x = d(L, ib, AW, LR); chk(f'aware-cap mean L={L} b#{ib}', x.mean(), mean); chk(f'aware-cap wins L={L} b#{ib}', (x > 0).sum(), wins, 0)
a2c = allpairs(A2, LR); a2a = allpairs(A2, AW)
chk('aware2-cap min', a2c.min(), 0.26); chk('aware2-cap max', a2c.max(), 5.09); chk('aware2-cap mean', a2c.mean(), 2.12)
chk('aware2-cap positive', (a2c > 0).sum(), 24, 0)
chk('aware2-aware min', a2a.min(), 0.26); chk('aware2-aware max', a2a.max(), 2.01); chk('aware2-aware mean', a2a.mean(), 0.81)
chk('aware2-aware positive', (a2a > 0).sum(), 24, 0)
print('all K=4 claims verified')
