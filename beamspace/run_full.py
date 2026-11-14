"""Full experiment grid (resumable).  usage: python run_full.py OUT.jsonl SECONDS DRAWS [scenarios] [Ks] [Ls]
Schemes per (scenario, K, L, draw):
  act_gG   : beamspace (lens) + active RIS, users split into groups of size G (time-shared), Ns = ceil(1.5 G) beams/slot
  pas_gG   : same with a passive RIS (no amplifier / amplifier noise)
  full_gG  : full-digital (all M = 64 beams, i.e. 64 RF chains) + active RIS
  *_all    : all K users served simultaneously
All designs use the proportional-fair phase objective (fair = 0.5) with gradient-refined RIS phases and RZF precoding."""
import os, sys, json, time, copy, itertools
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from multiprocessing import Pool
from bsmodel import *

SCN = {'S1': dict(dBR=50.0, dRU=20.0), 'S2': dict(dBR=20.0, dRU=10.0)}
PW = dict(eta_bs=0.4, P_rf=0.25, P_fix=1.0, eta_ris=0.4, P_el_act=0.02, P_el_pas=0.001)
FAIR = 0.5

def groups_by_correlation(f, g):
    K = f.shape[0]; fn = f / np.linalg.norm(f, axis=1, keepdims=True); Cc = np.abs(fn @ fn.conj().T)
    left = list(range(K)); out = []
    while left:
        grp = [left.pop(int(np.argmax(np.linalg.norm(f[left], axis=1))))]
        while len(grp) < g and left:
            j = int(np.argmin([Cc[i, grp].max() for i in left])); grp.append(left.pop(j))
        out.append(grp)
    return out

def configs(K, M=64):
    """(mode, group size, Ns) list. Ns > 24 is replaced by the full array (all M beams) to avoid the cost of greedy selection."""
    out = []
    for g in sorted({4, 8, K} & set(range(1, K + 1))):
        for Ns in sorted({min(M, int(np.ceil(m * g))) for m in (1.5, 3, 6)}):
            Ns = M if Ns > 24 else Ns
            out.append(('active', g, Ns))
        out.append(('active', g, M))
    out = sorted(set(out))
    pas = {('passive', 4 if K > 4 else K, 12), ('passive', K, 12 if K <= 8 else M)}
    return out + sorted(pas)

def serve(sc, g, mode, Ns_req):
    K = sc.K; A = np.sqrt(sc.Pris / sc.N); grps = groups_by_correlation(sc.f, min(g, K)); G = len(grps)
    Rk = np.zeros(K); Ps = []
    for grp in grps:
        sub = copy.copy(sc); sub.K = len(grp); sub.f = sc.f[grp]
        Ns = Ns_req
        S = np.arange(sc.M) if Ns >= sc.M else select_greedy(sub, Ns, np.inf, aware=True, cand=min(sc.M, max(24, 2 * Ns)))
        if mode == 'passive':
            sub.sv2 = 0.0; sub.Pris = np.inf
            d = design_refined(sub, S, np.inf, fair=FAIR, pgrid=np.array([1.0]), outer=2, steps=40); s, _ = evaluate(sub, d['S'], d['W'], d['ph'], 1.0, np.inf)
            P = sc.Pbs / PW['eta_bs'] + Ns * PW['P_rf'] + PW['P_fix'] + sc.N * PW['P_el_pas']
        else:
            d = design_refined(sub, S, A, fair=FAIR, outer=2, steps=40); s, Pr = evaluate(sub, d['S'], d['W'], d['ph'], d['p'], A)
            P = sc.Pbs / PW['eta_bs'] + Ns * PW['P_rf'] + PW['P_fix'] + Pr / PW['eta_ris'] + sc.N * PW['P_el_act']
        Rk[grp] = np.log2(1 + s); Ps.append(P)
    Rk = Rk / G; P = float(np.mean(Ps))
    return dict(sum=float(Rk.sum()), mn=float(Rk.min()), jain=float(Rk.sum() ** 2 / (K * (Rk ** 2).sum() + 1e-30)), P=P, EE=float(Rk.sum() / P), slots=G)

def task(a):
    t0 = time.time(); scn, K, L, dr = a
    sc = Scene(np.random.default_rng(100 + dr), K=K, L=L, **SCN[scn]); sc.Pris = 0.1
    out = {}
    for mode, g, Ns in configs(K):
        out[f"{mode[:3]}_g{g}_n{Ns}"] = serve(sc, g, mode, Ns)
    return dict(scn=scn, K=K, L=L, draw=dr, res=out, sec=round(time.time() - t0, 1))

if __name__ == '__main__':
    outf, budget, ndraw = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    scns = sys.argv[4].split(',') if len(sys.argv) > 4 else ['S2']
    Ks = [int(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else [4, 8, 16, 32]
    Ls = [int(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else [6, 16]
    done = set()
    if os.path.exists(outf):
        for ln in open(outf): r = json.loads(ln); done.add((r['scn'], r['K'], r['L'], r['draw']))
    todo = [(s, K, L, d) for d in range(ndraw) for K in Ks for s in scns for L in Ls if (s, K, L, d) not in done]
    t0 = time.time(); print('todo', len(todo), flush=True)
    with Pool(2) as pool, open(outf, 'a') as fo:
        it = pool.imap_unordered(task, todo)
        for r in it:
            fo.write(json.dumps(r) + '\n'); fo.flush()
            print('done', r['scn'], r['K'], r['L'], r['draw'], f'{time.time()-t0:.0f}s', flush=True)
            if time.time() - t0 > budget: pool.terminate(); break
