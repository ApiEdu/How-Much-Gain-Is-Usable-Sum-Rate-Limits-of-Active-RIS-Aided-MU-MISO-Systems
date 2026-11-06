"""Print the LaTeX rows of the K = 4 table (Table tab:k4 in the manuscript) from the result files."""
import analyze_k4 as A
import numpy as np

R = {L: A.load(L) for L in (100, 200)}
n = R[100].shape[0]
def cell(L, ib, s):
    x = R[L][:, ib, s]; return f'${x.mean():.2f}\\pm{A.se(x):.2f}$'
def dcell(L, ib, a, c):
    d = R[L][:, ib, a] - R[L][:, ib, c]
    return f'${d.mean():+.2f}$ ({int((d > 0).sum())}/{n})'.replace('$-', '$-').replace('-', '$-$') if False else f'${d.mean():+.2f}$ ({int((d > 0).sum())}/{n})'
rows = [('Grad., linear model', 3), ('Grad., linear-region cap', 1), ('Grad., compression-aware', 0),
        ('Grad., aware, start: cap$^\\ddagger$', 2)]
for name, s in rows:
    print(name + ' & ' + ' & '.join(cell(L, ib, s) for L in (100, 200) for ib in range(3)) + ' \\\\')
print('\\midrule')
for name, a, c in [('Grad. cap $-$ linear model', 1, 3), ('Compression-aware $-$ cap', 0, 1),
                   ('Aware, start: cap$^\\ddagger$ $-$ cap', 2, 1)]:
    print(name + ' & ' + ' & '.join(dcell(L, ib, a, c) for L in (100, 200) for ib in range(3)) + ' \\\\')
