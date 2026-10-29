# J98-style WMMSE/FP with the linear-region gain cap INSIDE the RIS update.
# Start: the 6000-iteration J98 solution with the cap applied afterwards (the "FP + gain cap" design).
# Each outer iteration: beamformer step as in fp_j98.py; RIS step solved by cvxpy with the per-element
# amplitude bound c_n = A*sqrt(rho_max/Pu_n), Pu_n = |G_n W|^2 + sv2 (W fixed inside the step).
import sys; sys.path.insert(0,'.')
import numpy as np, scipy.io as sio, cvxpy as cp, time
import fp_j98 as F
M,N,K=F.M,F.N,F.K; Pbs,Pris,sv2,s2,pmax=F.Pbs,F.Pris,F.sv2,F.s2,F.pmax; RHO=0.1
_x=cp.Variable(N,complex=True); _B=cp.Parameter((N,N),complex=True); _l=cp.Parameter(N,complex=True)
_c=cp.Parameter(N,nonneg=True); _U=cp.Parameter(N,nonneg=True)
_prob=cp.Problem(cp.Minimize(cp.sum_squares(_B@_x)+2*cp.real(_l@_x)),[cp.abs(_x)<=_c, cp.sum(cp.multiply(_U,cp.square(cp.abs(_x))))<=Pris])
def caps(G,W,A): return A*np.sqrt(RHO/(np.sum(np.abs(G@W)**2,1)+sv2))
def solve_cap(h,f,G,A,psi,W,iters):
    c=np.minimum(caps(G,W,A),pmax); psi=psi*np.minimum(1,c/np.maximum(np.abs(psi),1e-30))
    # post-hoc starting point must also satisfy the RIS power budget
    psi=F.project(G,psi,W)
    hist=[F.rate(h,f,G,psi,W)]; best=(hist[0],psi.copy(),W.copy())
    for it in range(iters):
        H=F.heq(h,f,G,psi); nz=s2+sv2*np.sum(np.abs(f*psi[None,:])**2,1)
        T=np.sum(np.abs(H@W)**2,1)+nz; u=np.diag(H@W)/T; e=1-np.abs(np.diag(H@W))**2/T; om=1/e
        Am=sum(om[k]*abs(u[k])**2*np.outer(H[k].conj(),H[k]) for k in range(K))
        v=np.stack([om[k]*H[k].conj()*u[k] for k in range(K)],1)
        D=G.conj().T@(np.abs(psi)[:,None]**2*G); Dh=(D+D.conj().T)/2
        w=cp.Variable((M,K),complex=True)
        obj=sum(cp.real(cp.quad_form(w[:,j],cp.Constant((Am+Am.conj().T)/2))) for j in range(K))-2*cp.real(cp.sum(cp.multiply(v.conj(),w)))
        cons=[cp.sum_squares(cp.abs(w))<=Pbs, sum(cp.real(cp.quad_form(w[:,j],cp.Constant(Dh))) for j in range(K))<=Pris-sv2*np.sum(np.abs(psi)**2)]
        try:
            cp.Problem(cp.Minimize(obj),cons).solve(solver=cp.CLARABEL)
            if w.value is not None: W=w.value
        except Exception: pass
        H=F.heq(h,f,G,psi); T=np.sum(np.abs(H@W)**2,1)+nz; u=np.diag(H@W)/T; e=1-np.abs(np.diag(H@W))**2/T; om=1/e
        Lam=np.zeros((N,N),complex); l=np.zeros(N,complex)
        for k in range(K):
            cc=h[k]@W; B=f[k][None,:]*(G@W).T
            Lam+=om[k]*abs(u[k])**2*(B.conj().T@B+sv2*np.diag(np.abs(f[k])**2))
            l+=om[k]*(abs(u[k])**2*(cc.conj()@B)-np.conj(u[k])*B[k])
        Lh=(Lam+Lam.conj().T)/2; lam,V=np.linalg.eigh(Lh); lam=np.maximum(lam,0)
        Bm=np.sqrt(lam)[:,None]*V.conj().T
        Ud=np.sum(np.abs(G@W)**2,1)+sv2; c=np.minimum(caps(G,W,A),pmax)
        _B.value=Bm; _l.value=l; _c.value=c; _U.value=Ud
        try:
            _prob.solve(solver=cp.CLARABEL)
            if _x.value is not None: psi=np.asarray(_x.value).ravel()
        except Exception: pass
        psi=psi*np.minimum(1,c/np.maximum(np.abs(psi),1e-30)); psi=F.project(G,psi,W)
        r=F.rate(h,f,G,psi,W); hist.append(r)
        if r>best[0]: best=(r,psi.copy(),W.copy())
    return best,hist
if __name__=="__main__":
    L=int(sys.argv[1]); d0,d1=int(sys.argv[2]),int(sys.argv[3]); ITERS=int(sys.argv[4]); bdb=[-5,0,5]
    for d in range(d0,d1+1):
        m=sio.loadmat(f"ch/ch_L{L}_d{d}.mat"); h,f,G=m['h'],m['f'],m['G']
        q=sio.loadmat(f"ch/fp6k_L{L}_d{d}.mat"); psi0=(q['p'].ravel()*np.exp(1j*q['theta'].ravel())); W0=q['w']
        for ib,b in enumerate(bdb):
            A=np.sqrt(10**(b/10)*Pris/N); t=time.time()
            (r,psi,W),hist=solve_cap(h,f,G,A,psi0.copy(),W0.copy(),ITERS)
            print(L,d,b,"linear-model rate start %.3f best %.3f last %.3f  time %.0fs"%(hist[0],r,hist[-1],time.time()-t),flush=True)
            sio.savemat(f"ch/fpcap{ITERS}_L{L}_d{d}_b{b+5}.mat",{'p':np.abs(psi),'theta':np.angle(psi),'w':W,'hist':np.array(hist)})
