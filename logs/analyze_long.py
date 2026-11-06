"""Long-budget check (run_long.m): K = 2, N = 256, 600 gradient iterations.

Reads logs/long_L<L>_d<draw>.log (lines  LONG <L> <draw> <beta_dB> <mode> <iteration> <rate> <J>)
and prints, for every (L, draw, beta), the symbol-level sum-rate of the capped design ("linreg") and of the
compression-aware design ("aware") after 150, 300 and 600 iterations, their difference (aware - linreg), and the
design-time objective J of the final iterate. Both designs start from their own closed-form design, as in
run_main150.m. Draws that are not present are skipped.
"""
import os, glob, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
CPS = (150, 300, 600)


def load():
    """Return {(L, draw, beta, mode, it): (rate, J)}."""
    data = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'long_L*_d*.log'))):
        for ln in open(f):
            if not ln.startswith('LONG '):
                continue
            t = ln.split()
            L, d, b, mode, it = int(t[1]), int(t[2]), int(t[3]), t[4], int(t[5])
            data[(L, d, b, mode, it)] = (float(t[6]), float(t[7]))
    return data


def main():
    data = load()
    keys = sorted({(L, d, b) for (L, d, b, m, it) in data})
    rows = collections.OrderedDict()
    for (L, d, b) in keys:
        print(f'\nL = {L} m, draw {d}, beta = {b:+d} dB   (symbol-level rate / design-time J)')
        for it in CPS:
            ra, ja = data[(L, d, b, 'aware', it)]
            rc, jc = data[(L, d, b, 'linreg', it)]
            print(f'  {it:4d} it   aware {ra:6.3f} / {ja:6.3f}   cap {rc:6.3f} / {jc:6.3f}   '
                  f'aware-cap: rate {ra - rc:+.3f}   J {ja - jc:+.3f}')
            rows[(L, d, b, it)] = dict(ra=ra, rc=rc, ja=ja, jc=jc)
        # change from 150 to 600 iterations
        ra0, ra1 = data[(L, d, b, 'aware', 150)][0], data[(L, d, b, 'aware', 600)][0]
        rc0, rc1 = data[(L, d, b, 'linreg', 150)][0], data[(L, d, b, 'linreg', 600)][0]
        print(f'  150 -> 600 it: aware {ra1 - ra0:+.3f}, cap {rc1 - rc0:+.3f}')
    return rows


def summary(rows=None):
    """Means over the available draws and counts used in the manuscript."""
    if rows is None:
        rows = main_quiet()
    out = {}
    for L in (100, 200):
        for b in (-5, 0):
            ds = sorted({d for (l, d, bb, it) in rows if l == L and bb == b})
            if not ds:
                continue
            for it in CPS:
                r = [rows[(L, d, b, it)] for d in ds]
                out[(L, b, it)] = dict(
                    n=len(ds),
                    cap=sum(x['rc'] for x in r) / len(r),
                    aware=sum(x['ra'] for x in r) / len(r),
                    dr=sum(x['ra'] - x['rc'] for x in r) / len(r),
                    dJ=sum(x['ja'] - x['jc'] for x in r) / len(r))
    return out


def main_quiet():
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return main()


if __name__ == '__main__':
    rows = main()
    print('\n=== means over draws (symbol-level rate; aware - cap) ===')
    for (L, b, it), v in summary(rows).items():
        print(f'L={L:3d} beta={b:+d} it={it:3d} n={v["n"]}  cap {v["cap"]:.3f}  aware {v["aware"]:.3f}  '
              f'aware-cap rate {v["dr"]:+.3f}  J {v["dJ"]:+.3f}')
    # J minus symbol-level rate of the final design ("design-time optimism"), 600 iterations
    print('\n=== J - symbol-level rate at 600 iterations (aware / cap) ===')
    for (L, d, b, it), x in rows.items():
        if it == 600:
            print(f'L={L} draw {d} beta={b:+d}: aware {x["ja"] - x["ra"]:+.3f}  cap {x["jc"] - x["rc"]:+.3f}  '
                  f'aware excess over cap {(x["ja"] - x["ra"]) - (x["jc"] - x["rc"]):+.3f}')
