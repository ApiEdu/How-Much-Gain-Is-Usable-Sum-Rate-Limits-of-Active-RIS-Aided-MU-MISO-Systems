"""Fill paper/main.tex.in with numbers from paper/numbers.json (+ validation.json, decouple.json), order the bibliography by
first citation, write paper/main.tex and compile it.   usage: python build_paper.py   (run make_figs.py first)"""
import json, re, subprocess, sys, os
import numpy as np

N = json.load(open('paper/numbers.json')); VAL = json.load(open('validation.json')); DEC = json.load(open('decouple.json'))
W = []   # warnings: claims in the template that the data no longer supports
def f(x, d=1): return f'{x:.{d}f}'
def g(x): return f'{x:g}'
V = {'n': str(N['n']), 'er6': g(N['effrank20_L6']), 'er16': g(N['effrank20_L16'])}
for L in (6, 16):
    for gg in (4, 8, 16):
        V[f'S_{L}_{gg}'] = f(N[f'K16_L{L}_g{gg}_sum'][0]); V[f'Mn_{L}_{gg}'] = f(N[f'K16_L{L}_g{gg}_mn'][0], 2); V[f'J_{L}_{gg}'] = f(N[f'K16_L{L}_g{gg}_jain'][0], 2)
        V[f'best_{L}_{gg}'] = f(N[f'K16_L{L}_g{gg}_bestsum'][0])
V['A16_all'] = V['S_16_16']; V['A16_g8'] = V['S_16_8']
V.update(EEr_mean=f(N['ee_ratio_all'][0], 2), EEr_min=f(N['ee_ratio_min'], 2), EEr_max=f(N['ee_ratio_max'], 2), npairs=str(N['n_pairs']),
         loss_mean=f(100 * N['rate_loss_mean']), loss_max=f(100 * N['rate_loss_max']),
         sens_lo=f(N['sens_min_mean'], 2), sens_hi=f(N['sens_max_mean'], 2), sens_min=f(N['sens_min_min'], 2),
         pas_K8L16=f(N['pas_sum_K8L16'][0]), act_K8L16=f(N['act_sum_K8L16'][0]), pas_max=f(N['pas_sum_max']))
K16 = N['K16_all']
V['ee16_16'] = f(N['ee_best']['16_16'][0], 2); V['ee16_6'] = f(N['ee_best']['16_6'][0], 2)
def cfg(k): m = re.match(r'act_g(\d+)_n(\d+)', k); return f'$g={m.group(1)}$, $N_s={m.group(2)}$'
V['ee16_16_cfg'] = cfg(N['ee_best']['16_16'][1]); V['ee16_6_cfg'] = cfg(N['ee_best']['16_6'][1])
V['ee64_g4_16'] = f(K16['16']['act_g4_n64']['EE'][0], 2); V['ee64_g4_6'] = f(K16['6']['act_g4_n64']['EE'][0], 2)
V['P12'] = f(K16['16']['act_g4_n12']['P'][0]); V['P64'] = f(K16['16']['act_g4_n64']['P'][0])
for L in (6, 16): V[f'gain64_{L}'] = f(100 * (K16[str(L)]['act_g4_n64']['sum'][0] / K16[str(L)]['act_g4_n12']['sum'][0] - 1))
fx = N['fixed_ratio']
V.update(fx_g4n6=f(fx['g4_n6'][0][0], 2), fx_g4n12=f(fx['g4_n12'][0][0], 2), fx_g4n12_min=f(fx['g4_n12'][1], 2), fx_g8n12=f(fx['g8_n12'][0][0], 2), fx_g8n12_min=f(fx['g8_n12'][1], 2))
V.update(val_n=str(VAL['n']), val_mean_abs=f(VAL['mean_abs'], 2), val_mean_rel_pct=f(100 * VAL['mean_rel'], 2), val_max_abs=f(VAL['max_abs'], 2), val_max_rel_pct=f(100 * VAL['max_rel'], 1), val_signed=f(VAL['mean_signed'], 2))
rho = [d['rho'] for d in DEC]; r5 = [d['rho'] for d in DEC if d['beta'] == -5]; rp = [d['rho'] for d in DEC if d['beta'] == 5]
V.update(dec_draws=str(len({d['draw'] for d in DEC})), dec_rho_min=f(min(rho), 3), dec_rho_p5=f(int(min(rp) * 1000) / 1000, 3), dec_same=str(sum(d['same_best'] for d in DEC)), dec_n=str(len(DEC)),
         dec_maxloss=f(max(d['loss'] for d in DEC), 2), dec_spread_lo=f(min(d['spread'] for d in DEC)), dec_spread_hi=f(max(d['spread'] for d in DEC)))
if min(rho) != min(r5): W.append('min rho is not at beta=-5 dB')
sb = N['sec_by_K']; V.update(sec4=f(sb['4'][0]), sec8=f(sb['8'][0]), sec16=f(sb['16'][0]), sec16max=f(sb['16'][1]))
# ---- checks of qualitative statements made in the template
if N['best_g_K16']['16'] != 8: W.append('best g for L=16 (N_s=1.5g) is not 8')
if N['best_g_K16']['6'] != 4: W.append('best g for L=6 (N_s=1.5g) is not 4')
for L, gg in ((16, 8), (6, 4)):
    bs = {x: N[f'K16_L{L}_g{x}_bestsum'][0] for x in (4, 8, 16)}
    if max(bs, key=bs.get) != gg: W.append(f'best-sum g for L={L} is {max(bs, key=bs.get)}, not {gg}')
for kk in ('16_16', '16_6'):
    if int(re.match(r'act_g\d+_n(\d+)', N['ee_best'][kk][1]).group(1)) > 12: W.append(f'best-EE config at {kk} uses more than 12 RF chains (template says 6-12)')
nonmono = any(K16[str(L)][f'act_g{gg}_n64']['sum'][0] < K16[str(L)][k]['sum'][0] for L in (6, 16) for gg, k in ((16, 'act_g16_n24'), (8, 'act_g8_n12'), (4, 'act_g4_n12')))
if not nonmono: W.append('sum-rate monotonic in N_s (template says it is not)')
if not (N['K16_L6_g4_jain'][0] < 0.9): W.append('Jain at L=6,g=4 >= 0.9 (template says poor fairness)')
if not (N['rate_loss_mean'] < 0.02): W.append('mean rate loss >= 2% (template says near zero)')

# ---- bibliography
BIB = {
 'long2021': 'R. Long, Y.-C. Liang, Y. Pei, and E. G. Larsson, ``Active reconfigurable intelligent surface-aided wireless communications,\'\' \\emph{IEEE Trans. Wireless Commun.}, vol. 20, no. 8, pp. 4962--4975, Aug. 2021.',
 'zhang2023': 'Z. Zhang, L. Dai, X. Chen, C. Liu, F. Yang, R. Schober, and H. V. Poor, ``Active RIS vs. passive RIS: Which will prevail in 6G?\'\' \\emph{IEEE Trans. Commun.}, vol. 71, no. 3, pp. 1707--1725, Mar. 2023.',
 'zeng2016': 'Y. Zeng and R. Zhang, ``Millimeter wave MIMO with lens antenna array: A new path division multiplexing paradigm,\'\' \\emph{IEEE Trans. Commun.}, vol. 64, no. 4, pp. 1557--1571, Apr. 2016.',
 'sayeed2013': 'A. Sayeed and J. Brady, ``Beamspace MIMO for high-dimensional multiuser communication at millimeter-wave frequencies,\'\' in \\emph{Proc. IEEE GLOBECOM}, Atlanta, GA, USA, Dec. 2013, pp. 3679--3684.',
 'wang2017': 'B. Wang, L. Dai, Z. Wang, N. Ge, and S. Zhou, ``Spectrum and energy efficient beamspace MIMO-NOMA for millimeter-wave communications using lens antenna array,\'\' \\emph{IEEE J. Sel. Areas Commun.}, vol. 35, no. 10, pp. 2370--2382, Oct. 2017.',
 'shen2019': 'W. Shen, X. Bu, X. Gao, C. Xing, and L. Hanzo, ``Beamspace precoding and beam selection for wideband millimeter-wave MIMO relying on lens antenna arrays,\'\' \\emph{IEEE Trans. Signal Process.}, vol. 67, no. 24, pp. 6301--6313, Dec. 2019.',
 'wanggc2020': 'Y. Wang, H. Lu, D. Zhao, and H. Sun, ``Energy efficiency optimization in IRS-enhanced mmWave systems with lens antenna array,\'\' in \\emph{Proc. IEEE GLOBECOM}, 2020, doi: 10.1109/GLOBECOM42002.2020.9348266.',
 'liu2021': 'P. Liu, Y. Li, W. Cheng, X. Gao, and X. Huang, ``Intelligent reflecting surface aided NOMA for millimeter-wave massive MIMO with lens antenna array,\'\' \\emph{IEEE Trans. Veh. Technol.}, 2021, doi: 10.1109/TVT.2021.3067938.',
 'pal2023': 'R. Pal and K. P. Sarawadekar, ``Distributed RIS-assisted mmWave multi-user MIMO beamspace system,\'\' \\emph{AEU -- Int. J. Electron. Commun.}, 2023, doi: 10.1016/j.aeue.2023.154560.',
 'abdullah2025': 'Z. Abdullah, M. R. Camana, A. B. M. Adam, C. K. Sheemar, E. Lagunas, and S. Chatzinotas, ``Swarm intelligence optimization of multi-RIS aided mmWave beamspace MIMO,\'\' in \\emph{Proc. IEEE VTC2025-Spring}, 2025, doi: 10.1109/VTC2025-Spring65109.2025.11174900.',
 'interdonato2024': 'G. Interdonato, F. Di Murro, C. D\'Andrea, G. Di Gennaro, and S. Buzzi, ``Approaching massive MIMO performance with reconfigurable intelligent surfaces: We do not need many antennas,\'\' arXiv:2203.07493 (v2), Oct. 2024.',
 'huang2024': 'A. Huang, X. Mu, L. Guo, and G. Zhu, ``Hybrid active-passive RIS transmitter enabled energy-efficient multi-user communications,\'\' \\emph{IEEE Trans. Wireless Commun.}, vol. 23, no. 9, pp. 10653--10666, Sep. 2024.',
 'yang2024': 'Y. Yang, Z. Lu, M. Li, R. Liu, and Q. Liu, ``A practical beamforming design for active RIS-assisted MU-MISO systems,\'\' arXiv:2401.03873, 2024.',
 'peng2024': 'Z. Peng, Z. Zhang, C. Pan, M. Di Renzo, O. A. Dobre, and J. Wang, ``Beamforming optimization for active RIS-aided multiuser communications with hardware impairments,\'\' arXiv:2402.10687, 2024.',
 'gavriilidis2025': 'P. Gavriilidis, D. Mishra, B. Smida, E. Basar, C. Yuen, and G. C. Alexandropoulos, ``Active reconfigurable intelligent surfaces: Circuit modeling and reflection amplification optimization,\'\' arXiv:2503.24093, 2025.',
 'jiao2022': 'R. Jiao, L. Dai, W. Wang, F. Lyu, N. Cheng, and X. Shen, ``Max-min fairness for beamspace MIMO-NOMA: From single-beam to multi-beam,\'\' \\emph{IEEE Trans. Wireless Commun.}, vol. 21, no. 2, pp. 739--752, Feb. 2022.',
 'rapp1991': 'C. Rapp, ``Effects of HPA-nonlinearity on a 4-DPSK/OFDM signal for a digital sound broadcasting system,\'\' in \\emph{Proc. 2nd Eur. Conf. Satell. Commun.}, Li\\`ege, Belgium, Oct. 1991, pp. 179--184.',
 'buss1952': 'J. J. Bussgang, ``Crosscorrelation functions of amplitude-distorted Gaussian signals,\'\' MIT Res. Lab. Electron., Cambridge, MA, USA, Tech. Rep. 216, Mar. 1952.',
 'kingma2015': 'D. P. Kingma and J. Ba, ``Adam: A method for stochastic optimization,\'\' in \\emph{Proc. Int. Conf. Learn. Represent. (ICLR)}, 2015, arXiv:1412.6980.',
 'peel2005': 'C. B. Peel, B. M. Hochwald, and A. L. Swindlehurst, ``A vector-perturbation technique for near-capacity multiantenna multiuser communication---Part I: Channel inversion and regularization,\'\' \\emph{IEEE Trans. Commun.}, vol. 53, no. 1, pp. 195--202, Jan. 2005.',
 'jain1984': 'R. Jain, D.-M. Chiu, and W. Hawe, ``A quantitative measure of fairness and discrimination for resource allocation in shared computer systems,\'\' Digital Equipment Corp., Tech. Rep. DEC-TR-301, Sep. 1984.',
 'koc': 'A. Koc and T. Le-Ngoc, ``Full-duplex mmWave massive MIMO systems: A joint hybrid precoding/combining and self-interference cancellation design,\'\' arXiv:2104.00191.',
 'qi2022': 'C. Qi, Q. Liu, X. Yu, and G. Y. Li, ``Hybrid precoding for mixture use of phase shifters and switches in mmWave massive MIMO,\'\' arXiv:2208.00714.',
 'hao2017': 'W. Hao, M. Zeng, Z. Chu, and S. Yang, ``Energy-efficient power allocation in millimeter wave massive MIMO with non-orthogonal multiple access,\'\' arXiv:1707.04520, 2017.',
}
tpl = open('paper/main.tex.in').read()
body = re.sub(r'@@(\w+)@@', lambda m: V[m.group(1)] if m.group(1) in V else m.group(0), tpl)
left = sorted(set(re.findall(r'@@(\w+)@@', body)) - {'BIB'})
assert not left, f'undefined placeholders: {left}'
order = []
for m in re.finditer(r'\\cite\{([^}]*)\}', body):
    for k in m.group(1).split(','):
        k = k.strip()
        if k not in order: order.append(k)
missing = [k for k in order if k not in BIB]; assert not missing, missing
unused = [k for k in BIB if k not in order]
bib = '\\begin{thebibliography}{99}\n' + '\n'.join(f'\\bibitem{{{k}}} {BIB[k]}' for k in order) + '\n\\end{thebibliography}'
body = body.replace('@@BIB@@', bib)
open('paper/main.tex', 'w').write(body)
print('placeholders filled; refs:', len(order), 'unused bib entries:', unused)
for w in W: print('WARNING (template claim not supported):', w)
if '--no-compile' not in sys.argv:
    for _ in range(2):
        r = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'main.tex'], cwd='paper', capture_output=True, text=True)
        if r.returncode: print(r.stdout[-3000:]); sys.exit(1)
    print('compiled paper/main.pdf')
