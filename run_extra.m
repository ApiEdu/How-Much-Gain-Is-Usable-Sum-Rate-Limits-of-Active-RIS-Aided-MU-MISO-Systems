% RUN_EXTRA  Gradient-optimised linear-region and compression-blind designs, same channels as run_study.
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000; giters=str2num(getenv('GITERS'));
Lset=str2num(getenv('LSET')); bdb=str2num(getenv('BDB')); ndraw=str2num(getenv('NDRAW')); tag=getenv('TAG');
for iL = 1:numel(Lset)
  L = Lset(iL);  up=[L,0,0; L+1,0,0];
  for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
  Rlr = nan(ndraw,numel(bdb));  Run = Rlr;  Raw = Rlr;
  for d = 1:ndraw
    rand('state',1000*L+d); randn('state',1000*L+d);
    [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
    crn = make_crn(K,N,S);
    ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax, ...
                  'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
    % compression-blind gradient design (linear model), independent of beta
    cf0 = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,1e3,pmax,crn);      % A irrelevant for .unaware
    gu = grad_design(ctx0, cf0.unaware, giters);
    for ib = 1:numel(bdb)
      A = sqrt(10^(bdb(ib)/10)*Pris/N);  ctx = ctx0; ctx.A = A;
      cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
      % (a) gradient design restricted to the linear region
      ctx.mode='linreg';
      d0 = cf.aware;  Pu = sum(abs(G*d0.w).^2,2)+sv2;  d0.p = max(min(d0.p, A*sqrt(rho_max./Pu)),1+1e-3);
      gl = grad_design(ctx, d0, giters);
      D = to_design(gl);  Pu = sum(abs(G*D.w).^2,2)+sv2;  D.p = min(D.p, A*sqrt(rho_max./Pu));
      sinr = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);  Rlr(d,ib) = sum(log2(1+sinr));
      % (b) compression-blind gradient design, re-scored on the compressed hardware
      ctx.mode='aware';  ind = gu;  ind = project_ris(ind, ctx);  D = to_design(ind);
      sinr = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);  Run(d,ib) = sum(log2(1+sinr));
      printf('L=%d draw=%d beta=%+d dB | Grad-linregion %.2f  Grad-unaware %.2f\n', L,d,bdb(ib),Rlr(d,ib),Run(d,ib)); fflush(stdout);
    end
    save(sprintf('res_%s_L%d.mat',tag,L),'Rlr','Run','bdb','L','d');
  end
end
