# WMMSE/FP alternating optimisation of (W, active-RIS psi) under the LINEAR active-RIS model (J98 / Zhang-Dai style).
import numpy as np, scipy.io as sio, cvxpy as cp, sys, glob
M,N,K=4,256,2; Pbs=9.9; Pris=0.1; sv2=1e-10; s2=1e-10; pmax=100.0
def heq(h,f,G,psi): return h+(f*psi[None,:])@G
def rate(h,f,G,psi,W):
    H=heq(h,f,G,psi); S=np.abs(H@W)**2
    nz=s2+sv2*np.sum(np.abs(f*psi[None,:])**2,1)
    sig=np.diag(S); itf=S.sum(1)-sig
    return np.sum(np.log2(1+sig/(itf+nz)))
def ris_power(G,psi,W): return np.sum(np.abs(psi[:,None]*(G@W))**2)+sv2*np.sum(np.abs(psi)**2)
def project(G,psi,W):
    psi=psi*np.minimum(1,pmax/np.maximum(np.abs(psi),1e-30))
    P=ris_power(G,psi,W)
    if P>Pris:  # shrink psi (power is quadratic in psi)
        psi=psi*np.sqrt(Pris/P)
    return psi
def solve(h,f,G,iters=1500):
    gh=h[0].conj()/np.linalg.norm(h[0]); ph=-np.angle(f[0]*(G@gh))+np.angle(h[0]@gh)
    psi=np.exp(1j*ph)*5.0
    H=heq(h,f,G,psi); W=np.linalg.pinv(H); W*=np.sqrt(Pbs/np.sum(np.abs(W)**2))
    psi=project(G,psi,W); best=(rate(h,f,G,psi,W),psi.copy(),W.copy()); hist=[best[0]]
    for it in range(iters):
        H=heq(h,f,G,psi); nz=s2+sv2*np.sum(np.abs(f*psi[None,:])**2,1)
        T=np.sum(np.abs(H@W)**2,1)+nz; u=np.diag(H@W)/T; e=1-np.abs(np.diag(H@W))**2/T; om=1/e
        A=sum(om[k]*abs(u[k])**2*np.outer(H[k].conj(),H[k]) for k in range(K))
        v=np.stack([om[k]*H[k].conj()*u[k] for k in range(K)],1)   # M x K
        D=G.conj().T@(np.abs(psi)[:,None]**2*G)
        w=cp.Variable((M,K),complex=True)
        Ah=np.linalg.cholesky(A+1e-30*np.eye(M)).conj().T if False else None
        obj=sum(cp.real(cp.quad_form(w[:,j],cp.Constant((A+A.conj().T)/2))) for j in range(K))-2*cp.real(cp.sum(cp.multiply(v.conj(),w)))
        Dh=(D+D.conj().T)/2
        cons=[cp.sum_squares(cp.abs(w))<=Pbs, sum(cp.real(cp.quad_form(w[:,j],cp.Constant(Dh))) for j in range(K))<=Pris-sv2*np.sum(np.abs(psi)**2)]
        try:
            cp.Problem(cp.Minimize(obj),cons).solve(solver=cp.CLARABEL)
            if w.value is not None: W=w.value
        except Exception as ex: pass
        # psi update
        H=heq(h,f,G,psi); T=np.sum(np.abs(H@W)**2,1)+nz; u=np.diag(H@W)/T; e=1-np.abs(np.diag(H@W))**2/T; om=1/e
        Lam=np.zeros((N,N),complex); l=np.zeros(N,complex)
        for k in range(K):
            c=h[k]@W                       # (K,) direct part c_kj
            B=f[k][None,:]*(G@W).T         # (K,N): beta_kj
            Lam+=om[k]*abs(u[k])**2*(B.conj().T@B+sv2*np.diag(np.abs(f[k])**2))
            l+=om[k]*(abs(u[k])**2*(c.conj()@B)-np.conj(u[k])*B[k])
        U=np.diag(np.sum(np.abs(G@W)**2,1)+sv2).astype(float)
        sol=lambda mu: np.linalg.solve(Lam+mu*U+1e-30*np.eye(N),-l.conj())
        x=sol(0.0)
        if np.real(x.conj()@U@x)>Pris:
            lo,hi=0.0,1.0
            while np.real(sol(hi).conj()@U@sol(hi))>Pris: lo,hi=hi,hi*10
            for _ in range(60):
                mid=(lo+hi)/2
                if np.real(sol(mid).conj()@U@sol(mid))>Pris: lo=mid
                else: hi=mid
            x=sol(hi)
        psi=project(G,x,W)
        r=rate(h,f,G,psi,W); hist.append(r)
        if r>best[0]: best=(r,psi.copy(),W.copy())
    return best,hist
if __name__=="__main__":
    L=int(sys.argv[1]); d0,d1=int(sys.argv[2]),int(sys.argv[3]); ITERS=int(sys.argv[4]) if len(sys.argv)>4 else 1500; PFX=sys.argv[5] if len(sys.argv)>5 else 'fp'
    for d in range(d0,d1+1):
        m=sio.loadmat(f"ch/ch_L{L}_d{d}.mat"); h,f,G=m['h'],m['f'],m['G']
        (r,psi,W),hist=solve(h,f,G,iters=ITERS)
        print(L,d,"linear-model rate %.2f (iter0 %.2f) Pris %.3f Pbs %.3f"%(r,hist[0],ris_power(G,psi,W)/Pris,np.sum(np.abs(W)**2)/Pbs),flush=True)
        sio.savemat(f"ch/{PFX}_L{L}_d{d}.mat",{'p':np.abs(psi),'theta':np.angle(psi),'w':W,'rate_lin':r})
