"""Beamspace (lens-antenna) MU system aided by an active RIS with Rapp-type amplifiers.

Model (downlink, direct BS-user link blocked):
  x  = U[:,S] Wb s                      BS: lens (unitary DFT U), only the beams in S have an RF chain
  u  = G x + v                          RIS input, v ~ CN(0, sv2 I)
  y  = ph * g(u)                        g = Rapp amplifier (small-signal gain p, output saturation amplitude A)
  r  = f y + n                          users, n ~ CN(0, s2 I)
Bussgang: g(u_n) = kappa_n u_n + d_n,  d_n uncorrelated with u_n, power sd2_n (treated as independent across n).
"""
import numpy as np
from numpy.polynomial.laguerre import laggauss

_T, _W = laggauss(96)          # Gauss-Laguerre nodes/weights for t ~ Exp(1)


def rapp_bussgang(Pu, p, A):
    """Bussgang gain kappa, distortion power sd2 and mean output power Pout of the Rapp (s=2) limiter
    g(u)=p u/(1+(p|u|/A)^4)^(1/4) for u ~ CN(0,Pu). Vectorised over elements."""
    Pu = np.asarray(Pu, float); p = np.broadcast_to(np.asarray(p, float), Pu.shape)
    if np.isinf(A):
        return p.astype(float), np.zeros_like(Pu), p**2 * Pu
    rho = p**2 * Pu / A**2
    r = np.outer(rho, _T)
    kI = ((1 + r**2) ** -0.25 * _T) @ _W
    oI = ((1 + r**2) ** -0.5 * _T) @ _W
    kappa = p * kI
    Pout = p**2 * Pu * oI
    sd2 = np.maximum(Pout - kappa**2 * Pu, 0)
    return kappa, sd2, Pout


def steer(n, sin_t):
    return np.exp(1j * np.pi * np.arange(n)[:, None] * np.atleast_1d(sin_t)[None, :])


def dft_lens(M):
    """Unitary DFT of a lens array: column m is the beam towards sin(phi)=(2m-(M-1))/M."""
    m = np.arange(M); psi = (2 * m - (M - 1)) / M
    return steer(M, psi) / np.sqrt(M)


class Scene:
    def __init__(self, rng, M=64, N=128, K=4, L=6, fc=28e9, dBR=50.0, dRU=20.0, Pbs=1.0, Pris=1e-2,
                 sv2=2e-12, s2=2e-12, main_frac=0.7, min_sep=None):
        lam = 3e8 / fc
        pl = lambda d: (lam / (4 * np.pi * d)) ** 2
        self.M, self.N, self.K, self.L = M, N, K, L
        self.Pbs, self.Pris, self.sv2, self.s2 = Pbs, Pris, sv2, s2
        # BS->RIS: L paths, path 0 carries main_frac of the power, angles at both ends random
        pw = np.r_[main_frac, (1 - main_frac) * rng.dirichlet(np.ones(L - 1))]
        c = np.sqrt(pw) * np.exp(2j * np.pi * rng.random(L))
        phi = rng.uniform(-0.85, 0.85, L); psi = rng.uniform(-0.85, 0.85, L)    # sin(angle)
        self.phi = phi
        self.G = np.sqrt(pl(dBR) * N * M) * (steer(N, psi) * c) @ steer(M, phi).conj().T / np.sqrt(N * M)
        # RIS->users: LoS, well separated directions
        if min_sep is None: min_sep = min(0.12, 0.6 * 1.7 / K)
        while True:
            chi = np.sort(rng.uniform(-0.85, 0.85, K))
            if K == 1 or np.min(np.diff(chi)) > min_sep: break
        g = np.exp(2j * np.pi * rng.random(K))
        self.f = np.sqrt(pl(dRU)) * (steer(N, chi) * g).T.conj()
        self.U = dft_lens(M)
        self.Gb = self.G @ self.U            # N x M beamspace BS->RIS channel


def zf(Heq, Pbs, alpha=0.0):
    """ZF (alpha=0) or regularised ZF precoder; alpha is relative to the mean eigenvalue of Heq Heq^H."""
    H2 = Heq @ Heq.conj().T
    H2 = H2 + max(alpha, 1e-9) * np.trace(H2).real / H2.shape[0] * np.eye(H2.shape[0])
    W = Heq.conj().T @ np.linalg.inv(H2)
    return W * np.sqrt(Pbs) / np.linalg.norm(W)


def phases(sc, S, iters=6, p=1.0):
    """RIS phases by alternating maximisation of sum_k Re(f_k diag(ph) Gb w_k) with ZF precoder."""
    Gb = sc.Gb[:, S]
    ph = np.exp(-1j * np.angle(sc.f[0]) - 1j * np.angle(Gb.sum(1)))
    for _ in range(iters):
        W = zf((sc.f * ph) @ Gb, sc.Pbs)
        z = (sc.f.T * (Gb @ W)).sum(1)          # sum_k f_kn (Gb w_k)_n
        ph = np.exp(-1j * np.angle(z))
    return ph


def evaluate(sc, S, Wb, ph, p, A):
    """Analytic (Bussgang) per-user SINR, RIS output power."""
    Gb = sc.Gb[:, S]
    X = Gb @ Wb
    Pu = (np.abs(X) ** 2).sum(1) + sc.sv2
    kap, sd2, Pout = rapp_bussgang(Pu, p, A)
    Heff = (sc.f * (kap * ph)) @ Gb
    T = Heff @ Wb
    sig = np.abs(np.diag(T)) ** 2
    intf = (np.abs(T) ** 2).sum(1) - sig
    noise = (np.abs(sc.f) ** 2) @ (kap**2 * sc.sv2 + sd2) + sc.s2
    return sig / (intf + noise), Pout.sum()


def design(sc, S, A, aware=True, pgrid=None, ph=None, iters=2, alphas=(0.0, 0.1, 0.3, 1.0)):
    """Gain p (uniform) and ZF precoder for beam set S.
    aware=True : p chosen on the compression-aware rate subject to the true RIS power budget.
    aware=False: p chosen by the LINEAR power budget (what a compression-blind design does)."""
    S = np.asarray(S)
    if ph is None: ph = phases(sc, S)
    Gb = sc.Gb[:, S]
    if pgrid is None: pgrid = np.logspace(0, 2, 16)
    def build(p, Ae, al=0.0):
        W = zf((sc.f * (p * ph)) @ Gb, sc.Pbs, al)
        for _ in range(iters):
            Pu = (np.abs(Gb @ W) ** 2).sum(1) + sc.sv2
            kap, _, _ = rapp_bussgang(Pu, p, Ae)
            W = zf((sc.f * (kap * ph)) @ Gb, sc.Pbs, al)
        return W
    if not aware:
        best = None
        for p in pgrid:                       # largest p whose LINEAR RIS output power fits the budget
            W = build(p, np.inf)
            _, Pr = evaluate(sc, S, W, ph, p, np.inf)
            if Pr <= sc.Pris: best = (p, W)
        if best is None: best = (1.0, build(1.0, np.inf))
        p, W = best
        return dict(S=S, ph=ph, p=p, W=W)
    best, br = None, -1
    for al in alphas:
        for p in pgrid:
            W = build(p, A, al)
            s, Pr = evaluate(sc, S, W, ph, p, A)
            r = np.log2(1 + s).sum()
            if Pr <= sc.Pris and r > br: best, br = (p, W), r
    if best is None: best = (1.0, build(1.0, A))
    p, W = best
    return dict(S=S, ph=ph, p=p, W=W)


def rate_analytic(sc, d, A):
    s, _ = evaluate(sc, d['S'], d['W'], d['ph'], d['p'], A)
    return np.log2(1 + s).sum()


def rate_symbol(sc, d, A, T=40000, seed=0):
    """Symbol-level ground truth: Gaussian symbols through the actual Rapp amplifier."""
    rng = np.random.default_rng(seed)
    K, N = sc.K, sc.N
    cn = lambda *s: (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2)
    S = d['S']; Ub = sc.U[:, S]
    s = cn(K, T); x = Ub @ (d['W'] @ s)
    u = sc.G @ x + np.sqrt(sc.sv2) * cn(N, T)
    p = d['p']
    y = (p * u) / (1 + (p * np.abs(u) / A) ** 4) ** 0.25 * d['ph'][:, None] if np.isfinite(A) else p * u * d['ph'][:, None]
    r = sc.f @ y + np.sqrt(sc.s2) * cn(K, T)
    a = (r * s.conj()).mean(1)
    res = r - a[:, None] * s
    sinr = np.abs(a) ** 2 / (np.abs(res) ** 2).mean(1)
    return np.log2(1 + sinr).sum()


def cascaded_energy(sc, ph):
    return (np.abs((sc.f * ph) @ sc.Gb) ** 2).sum(0)


def select_energy(sc, Ns, ph):
    return np.sort(np.argsort(-cascaded_energy(sc, ph))[:Ns])


def select_greedy(sc, Ns, A, aware=True, cand=None, pgrid=None):
    """Forward greedy beam selection on the aware (or linear-model) rate; starts from the K strongest beams."""
    if cand is None: cand = min(sc.M, max(24, 2 * Ns))
    ph0 = phases(sc, np.arange(sc.M))
    en = cascaded_energy(sc, ph0)
    order = np.argsort(-en)
    S = list(order[:sc.K]); pool = list(order[sc.K:cand])
    pg = np.logspace(0, 2, 8) if pgrid is None else pgrid
    def score(S_):
        d = design(sc, S_, A, aware=True, pgrid=pg)
        return rate_analytic(sc, d, A) if aware else rate_analytic(sc, d, np.inf)
    while len(S) < Ns and pool:
        sc_ = [score(S + [b]) for b in pool]
        j = int(np.argmax(sc_)); S.append(pool.pop(j))
    return np.sort(S)
