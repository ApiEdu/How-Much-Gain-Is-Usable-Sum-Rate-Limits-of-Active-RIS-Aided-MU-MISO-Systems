#!/usr/bin/env python3
"""Summarise the release-the-cap experiment (warm_L*.log) and the rho_max / p_max sweeps (rho_*.log, pmax*_*.log).
Cross-checks stage 1 against the logged 50-iteration capped designs (main150_L*.log) and prints LaTeX rows."""
import re, glob, math, collections, sys
sys.argv = ['x']
import analyze_main150 as A
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
def mean_se(x):
    n = len(x)
    m = sum(x) / n
    if n < 2:
        return m, float('nan')
    v = sum((a - m) ** 2 for a in x) / (n - 1)
    return m, math.sqrt(v / n)

# ------------------------------------------------------------ release test
W = collections.defaultdict(dict)   # (L,b) -> {(tag,it): {d: val}}
for L in (100, 200):
    try:
        for ln in open(S + 'warm_L%d.log' % L):
            m = re.match(r'WARM (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
            if m:
                l, d, b, tag, it, v = int(m[1]), int(m[2]), int(m[3]), m[4], int(m[5]), float(m[6])
                W[(l, b)].setdefault((tag, it), {})[d] = v
    except FileNotFoundError:
        pass
# logged references (draws 1-4)
M = collections.defaultdict(dict)
for L in (100, 200):
    for ln in open(S + 'main150_L%d.log' % L):
        m = re.match(r'M150 (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
        if m:
            l, d, b, mode, it, v = int(m[1]), int(m[2]), int(m[3]), m[4], int(m[5]), float(m[6])
            M[(l, b, mode, it)][d] = v
print('== stage-1 reproduction (cap after 50 it vs logged linreg 50), per draw')
mx = 0
for (l, b), dct in sorted(W.items()):
    for d, v in sorted(dct.get(('cap', 50), {}).items()):
        ref = A.val('linreg', l, d, b, 50)
        if ref == ref:
            print(' L%d b%+d d%d new %.3f logged %.3f diff %+.3f' % (l, b, d, v, ref, v - ref))
            mx = max(mx, abs(v - ref))
print(' max |diff| =', mx)
print('\n== release test: totals of 150 iterations = 50 (capped start) + 100')
print(' L   beta  n | cap kept | cap released | released-kept (pos/n) | cold aware150 | cold cap150 | released-coldcap')
rows = []
for (l, b), dct in sorted(W.items()):
    ck = dct.get(('capcont', 100), {})
    aw = dct.get(('awarewarm', 100), {})
    ds = sorted(set(ck) & set(aw))
    if not ds:
        continue
    x = [ck[d] for d in ds]
    y = [aw[d] for d in ds]
    diff = [aw[d] - ck[d] for d in ds]
    cold_a = [A.val('aware', l, d, b, 150) for d in ds]
    cold_c = [A.val('linreg', l, d, b, 150) for d in ds]
    mx_, sx = mean_se(x); my, sy = mean_se(y); md, sd = mean_se(diff)
    pos = sum(1 for a in diff if a > 0)
    ca = mean_se(cold_a)[0] if cold_a else float('nan')
    cc = mean_se(cold_c)[0] if cold_c else float('nan')
    print(' %d  %+d  %d | %.2f+-%.2f | %.2f+-%.2f | %+.2f+-%.2f (%d/%d) | %.2f | %.2f | %+.2f' %
          (l, b, len(ds), mx_, sx, my, sy, md, sd, pos, len(ds), ca, cc, my - cc))
    rows.append((l, b, len(ds), mx_, sx, my, sy, md, sd, pos, ca, cc, diff))
    print('    per-draw released-kept:', ' '.join('%+.2f' % a for a in diff))
    h = [aw[d] - dct.get(('awarewarm', 50), {}).get(d, float('nan')) for d in ds]

# gain statistics
print('\n== gain statistics (PST lines): mean p, max p, fraction of elements at >99% of p_max')
for L in (100, 200):
    try:
        agg = collections.defaultdict(list)
        for ln in open(S + 'warm_L%d.log' % L):
            m = re.match(r'PST (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)', ln)
            if m:
                agg[(int(m[3]), m[4], int(m[5]))].append((float(m[6]), float(m[7]), float(m[8])))
        for k, v in sorted(agg.items()):
            print(' L%d beta %+d %-9s it %3d : mean p %.1f max p %.1f frac>=0.99pmax %.3f (n=%d)' %
                  (L, k[0], k[1], k[2], sum(a[0] for a in v) / len(v), max(a[1] for a in v), sum(a[2] for a in v) / len(v), len(v)))
    except FileNotFoundError:
        pass

# ------------------------------------------------------------ rho / pmax sweeps
print('\n== rho_max / p_max sweeps (cap design, 150 it, beta=0)')
R = collections.defaultdict(dict)
files = glob.glob(S + 'rho_L*.log') + glob.glob(S + 'pmax*_L*.log')
for f in files:
    for ln in open(f):
        m = re.match(r'RHO (\d+) (\d+) (-?\d+) ([\d.]+) (\d+) (\d+) ([\d.]+)', ln)
        if m:
            l, d, b, rho, pm, it, v = int(m[1]), int(m[2]), int(m[3]), float(m[4]), int(m[5]), int(m[6]), float(m[7])
            R[(l, b, rho, pm)][d] = v
base = {}
for L in (100, 200):
    for d in (1, 2, 3):
        v = A.val('linreg', L, d, 0, 150)
        R[(L, 0, 0.1, 100)][d] = v           # rho 0.1 / pmax 100 from the main run (same code, same seeds)
        base[(L, d)] = v
for l in (100, 200):
    for key in sorted(k for k in R if k[0] == l):
        dct = R[key]
        ds = sorted(d for d in dct if d <= 3)
        if not ds:
            continue
        x = [dct[d] for d in ds]
        m_, s_ = mean_se(x)
        print(' L%d beta %+d rho_max %.2f p_max %3d : %.2f +- %.2f (n=%d)  per-draw %s' % (l, key[1], key[2], key[3], m_, s_, len(ds), ' '.join('%.2f' % a for a in x)))
