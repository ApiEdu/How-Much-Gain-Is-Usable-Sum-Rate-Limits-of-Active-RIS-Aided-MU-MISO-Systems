#!/usr/bin/env python3
"""Design-time objective (J160, on the S=160 design samples), model rate on 2000 fresh samples (Rf), and symbol-level score
for the release-the-cap test (warm_L*.log). Produces per-condition means and counts; dumps ext_dj.json."""
import re, json, math, collections
import os
S = os.path.dirname(os.path.abspath(__file__)) + '/'   # logs and scripts live in the same folder (logs/)
DJ, SC = {}, {}
for L in (100, 200):
    for ln in open(S + 'warm_L%d.log' % L):
        m = re.match(r'DJ (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+) ([\d.]+) ([\d.]+)', ln)
        if m: DJ[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = (float(m[6]), float(m[7]), float(m[8]))
        m = re.match(r'WARM (\d+) (\d+) (-?\d+) (\w+) (\d+) ([\d.]+)', ln)
        if m: SC[(int(m[1]), int(m[3]), int(m[2]), m[4], int(m[5]))] = float(m[6])
def ms(x):
    n = len(x); m = sum(x) / n
    return m, (math.sqrt(sum((a - m) ** 2 for a in x) / (n - 1) / n) if n > 1 else float('nan'))
out = {}
print('L beta | dJ160 (rel-kept) | dRfresh (rel-kept) | dScore (rel-kept) | gap J160-Rfresh: cap-kept, released | #designs J160 rel>kept | start->rel50 J160 unchanged')
tot = collections.Counter()
for L in (100, 200):
    for b in (-5, 0, 5):
        ds = [d for d in (1, 2, 3, 4) if (L, b, d, 'awarewarm', 100) in DJ and (L, b, d, 'capcont', 100) in DJ]
        dJ = [DJ[(L, b, d, 'awarewarm', 100)][0] - DJ[(L, b, d, 'capcont', 100)][0] for d in ds]
        dR = [DJ[(L, b, d, 'awarewarm', 100)][1] - DJ[(L, b, d, 'capcont', 100)][1] for d in ds]
        dS = [SC[(L, b, d, 'awarewarm', 100)] - SC[(L, b, d, 'capcont', 100)] for d in ds]
        gk = [DJ[(L, b, d, 'capcont', 100)][0] - DJ[(L, b, d, 'capcont', 100)][1] for d in ds]
        ga = [DJ[(L, b, d, 'awarewarm', 100)][0] - DJ[(L, b, d, 'awarewarm', 100)][1] for d in ds]
        unch = sum(1 for d in ds if DJ[(L, b, d, 'awarewarm', 50)][0] == DJ[(L, b, d, 'cap', 50)][0])
        posJ = sum(1 for x in dJ if x > 0)
        out['%d_%d' % (L, b)] = dict(n=len(ds), dJ=ms(dJ), dR=ms(dR), dS=ms(dS), gap_kept=ms(gk), gap_rel=ms(ga), posJ=posJ, unch50=unch,
                                    posR=sum(1 for x in dR if x > 0), posS=sum(1 for x in dS if x > 0))
        print('%d %+d n=%d | %+.2f | %+.2f | %+.2f | %+.2f, %+.2f | %d/%d | %d/%d' % (L, b, len(ds), ms(dJ)[0], ms(dR)[0], ms(dS)[0], ms(gk)[0], ms(ga)[0], posJ, len(ds), unch, len(ds)))
        tot['n'] += len(ds); tot['posJ'] += posJ; tot['unch'] += unch
        tot['posR'] += out['%d_%d' % (L, b)]['posR']; tot['posS'] += out['%d_%d' % (L, b)]['posS']
# split: designs where score(rel) > score(kept)
print(dict(tot))
# overfit gap for all cap-start designs
g0 = [DJ[k][0] - DJ[k][1] for k in DJ if k[3] == 'cap']
print('cap-start gap J160 - Rfresh over', len(g0), 'designs: mean %.3f min %.3f max %.3f' % (sum(g0) / len(g0), min(g0), max(g0)))
# Bussgang-model fresh rate vs symbol score, awarewarm100 / capcont100
for tag in ('capcont', 'awarewarm'):
    e = [DJ[k][1] - SC[k] for k in DJ if k[3] == tag and k[4] == 100]
    print(tag, 'model(fresh) - symbol score at 100 its: mean %.3f min %.3f max %.3f n=%d' % (sum(e) / len(e), min(e), max(e), len(e)))
json.dump(out, open(S + 'ext_dj.json', 'w'), indent=1)
