#!/usr/bin/env python3
"""Recompute from the logs every number quoted in the new subsections (release test, rho_max / p_max sweeps) and assert it
matches the rounded value written in the manuscript. Prints OK per claim; raises on mismatch."""
import re, glob, sys
sys.argv = ['x']
import analyze_main150 as A
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
DJ, SC, PST = {}, {}, {}
for L in (100, 200):
    for ln in open(S + 'warm_L%d.log' % L):
        m = re.match(r'DJ (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)', ln)
        if m: DJ[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = tuple(float(m[i]) for i in (6, 7, 8))
        m = re.match(r'WARM (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
        if m: SC[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = float(m[6])
        m = re.match(r'PST (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)', ln)
        if m: PST[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = tuple(float(m[i]) for i in (6, 7, 8))
keys = [(L, b, d) for L in (100, 200) for b in (-5, 0, 5) for d in (1, 2, 3, 4)]
def chk(name, val, expect, tol=0.005):
    ok = abs(val - expect) <= tol
    print(('OK  ' if ok else 'FAIL'), name, round(val, 4), 'expected', expect)
    assert ok, name
ds = [SC[k + ('awarewarm', 100)] - SC[k + ('capcont', 100)] for k in keys]
chk('released-kept min', min(ds), -0.49); chk('released-kept max', max(ds), 1.06); chk('released-kept mean', sum(ds) / 24, 0.10, 0.005)
chk('positive total', sum(x > 0 for x in ds), 10, 0); chk('positive L100', sum(x > 0 for x in ds[:12]), 8, 0); chk('positive L200', sum(x > 0 for x in ds[12:]), 2, 0)
chk('J lower at 200m', sum(DJ[k + ('awarewarm', 100)][0] < DJ[k + ('capcont', 100)][0] for k in keys if k[0] == 200), 7, 0)
chk('J higher at 100m', sum(DJ[k + ('awarewarm', 100)][0] > DJ[k + ('capcont', 100)][0] for k in keys if k[0] == 100), 12, 0)
chk('no improvement in 50 its', sum(DJ[k + ('awarewarm', 50)][0] == DJ[k + ('cap', 50)][0] for k in keys), 4, 0)
chk('kept improved all', sum(DJ[k + ('capcont', 50)][0] > DJ[k + ('cap', 50)][0] for k in keys), 24, 0)
for L in (100, 200):
    for b in (-5, 0, 5):
        ref = [A.val('linreg', L, d, b, 50) for d in (1, 2, 3, 4)]
        assert max(abs(SC[(L, b, d, 'cap', 50)] - r) for d, r in zip((1, 2, 3, 4), ref)) < 0.0006
print('OK   stage-1 reproduction exact for 24 designs')
def mean_gap(L, tag, f): 
    v = [DJ[(L, b, d, tag, 100)][0] - DJ[(L, b, d, tag, 100)][1] for b in (-5, 0, 5) for d in (1, 2, 3, 4)]
    return v
rel100 = [sum(DJ[(100, b, d, 'awarewarm', 100)][0] - DJ[(100, b, d, 'awarewarm', 100)][1] for d in (1, 2, 3, 4)) / 4 for b in (-5, 0, 5)]
kep100 = [sum(DJ[(100, b, d, 'capcont', 100)][0] - DJ[(100, b, d, 'capcont', 100)][1] for d in (1, 2, 3, 4)) / 4 for b in (-5, 0, 5)]
chk('released gap min', min(rel100), 0.22, 0.006); chk('released gap max', max(rel100), 0.30, 0.006); chk('kept gap max', max(kep100), 0.06, 0.006)
mx = max(PST[k + ('awarewarm', 100)][1] for k in keys); chk('largest released gain', mx, 99.8, 0.06)
chk('max frac within 1% of pmax (released)', max(PST[k + ('awarewarm', 100)][2] for k in keys), 0.012, 0.0006)
# sweeps
R = {}
for f in glob.glob(S + 'rho_L*.log') + glob.glob(S + 'pmax*_L*.log'):
    for ln in open(f):
        m = re.match(r'RHO (\d+) (\d+) (-?\d+) ([\d.]+) (\d+) (\d+) ([\d.]+)', ln)
        if m: R[(int(m[1]), float(m[4]), int(m[5]), int(m[2]))] = float(m[7])
for L in (100, 200):
    for d in (1, 2, 3): R[(L, 0.1, 100, d)] = A.val('linreg', L, d, 0, 150)
assert len(R) == 2 * (6 * 3), len(R)
D = (1, 2, 3)
r = lambda L, rho, pm, d: R[(L, rho, pm, d)]
chk('rho .03 loss L100', -sum(r(100, .03, 100, d) - r(100, .1, 100, d) for d in D) / 3, 3.0, 0.05)
chk('rho .03 loss L200', -sum(r(200, .03, 100, d) - r(200, .1, 100, d) for d in D) / 3, 1.7, 0.05)
d03 = [r(L, .3, 100, d) - r(L, .1, 100, d) for L in (100, 200) for d in D]
chk('rho .3 lower in all', sum(x < 0 for x in d03), 6, 0); chk('rho .3 min drop', -max(d03), 0.32, 0.006); chk('rho .3 max drop', -min(d03), 2.13, 0.006)
d1 = [r(L, 1.0, 100, d) - r(L, .1, 100, d) for L in (100, 200) for d in D]
chk('rho 1 lower in 5', sum(x < 0 for x in d1), 5, 0); chk('rho 1 smallest drop', min(-x for x in d1 if x < 0), 0.05, 0.006); chk('rho 1 largest drop', max(-x for x in d1), 0.48, 0.006); chk('rho 1 gain', max(d1), 0.29, 0.006)
best = [max((0.03, .1, .3, 1.0), key=lambda x: r(L, x, 100, d)) for L in (100, 200) for d in D]
chk('rho .1 best in', sum(b == .1 for b in best), 5, 0)
p300 = [r(L, .1, 300, d) - r(L, .1, 100, d) for L in (100, 200) for d in D]; chk('pmax300 max abs', max(abs(x) for x in p300), 0.05, 0.006)
p30_100 = [r(100, .1, 30, d) - r(100, .1, 100, d) for d in D]; chk('pmax30 L100 min', -max(p30_100), 0.19, 0.006); chk('pmax30 L100 max', -min(p30_100), 0.28, 0.006)
p30_200 = [r(200, .1, 30, d) - r(200, .1, 100, d) for d in D]; chk('pmax30 L200 low', min(p30_200), 0.01, 0.006); chk('pmax30 L200 high', max(p30_200), 0.13, 0.006)
for L in (100, 200):
    diffs = [r(L, rho, 100, d) - A.val('unaware', L, d, 0, 150) for rho in (.03, .1, .3, 1.0) for d in D]
    print(L, 'cap-lin min %.2f max %.2f' % (min(diffs), max(diffs)))
chk('cap-lin L100 min', min(r(100, rho, 100, d) - A.val('unaware', 100, d, 0, 150) for rho in (.03, .1, .3, 1.0) for d in D), 3.1, 0.05)
chk('cap-lin L100 max', max(r(100, rho, 100, d) - A.val('unaware', 100, d, 0, 150) for rho in (.03, .1, .3, 1.0) for d in D), 6.6, 0.05)
chk('cap-lin L200 min', min(r(200, rho, 100, d) - A.val('unaware', 200, d, 0, 150) for rho in (.03, .1, .3, 1.0) for d in D), 9.7, 0.06)
chk('cap-lin L200 max', max(r(200, rho, 100, d) - A.val('unaware', 200, d, 0, 150) for rho in (.03, .1, .3, 1.0) for d in D), 11.9, 0.05)
for L, exp in ((100, 0.30), (200, -0.20)):
    chk('aware-cap on draws 1-3 L%d' % L, sum(A.val('aware', L, d, 0, 150) - r(L, .1, 100, d) for d in D) / 3, exp, 0.006)
print('ALL CLAIMS VERIFIED')
