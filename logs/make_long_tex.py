"""Print the rows of the long-budget table of the manuscript (tab:long) from logs/long_L*_d*.log.

Columns: condition | cap 150 | cap 600 | aware 150 | aware 600 | aware-cap at 150, 300, 600 (mean over draws)
         | aware-cap at 600 per draw (draw 1 / draw 2) | J difference at 600 (mean).
Run without -I:  python3 logs/make_long_tex.py
"""
import analyze_long as A


def sgn(x, nd=2):
    """Signed number in math mode; rounds half away from zero on the 3-decimal value (logs have 3 decimals)."""
    from decimal import Decimal, ROUND_HALF_UP
    q = Decimal(f'{abs(x):.3f}').quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)
    neg = x < 0 and q != 0
    return ('$-' if neg else '$+') + str(q) + '$'


def lev(x, nd=2):
    """Unsigned level, rounded half up on the 3-decimal value (same rule as verify_long_claims.py)."""
    from decimal import Decimal, ROUND_HALF_UP
    return str(Decimal(f'{x:.3f}').quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def main():
    rows = A.main_quiet()
    S = A.summary(rows)
    out = []
    for L in (100, 200):
        for b in (-5, 0):
            c150, c600 = S[(L, b, 150)], S[(L, b, 600)]
            d = [rows[(L, dr, b, 600)] for dr in (1, 2)]
            per = ' / '.join(sgn(x['ra'] - x['rc']) for x in d)
            lab = f'{L}~m, $\\beta={b:+d}$~dB'.replace('+0', '0')
            out.append(' & '.join([
                lab, lev(c150["cap"]), lev(c600["cap"]), lev(c150["aware"]), lev(c600["aware"]),
                sgn(S[(L, b, 150)]['dr']), sgn(S[(L, b, 300)]['dr']), sgn(S[(L, b, 600)]['dr']), per,
                sgn(c600['dJ'])]) + r' \\')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
