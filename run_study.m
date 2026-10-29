% RUN_STUDY  Compression-aware active-RIS design study (MU-MISO, Rapp amplifiers).
% Environment: LSET (user positions, m), BDB (saturation-to-budget ratios, dB), NDRAW, TAG.
% beta = A_sat^2 / (P_RIS_max/N): amplifier saturation power relative to the average per-element budget.
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
popsize=50; gens=str2num(getenv("GENS")); tsize=5; mrate=0.1; pmax=100; rho_max=0.1; S=160; Tsim=40000; giters=str2num(getenv("GITERS"));
Lset=str2num(getenv('LSET')); bdb=str2num(getenv('BDB')); ndraw=str2num(getenv('NDRAW')); tag=getenv('TAG');
schemes = {'CF-unaware','CF-aware','GA-unaware','GA-linregion','GA-aware','Grad-aware'};
nS = numel(schemes);
for iL = 1:numel(Lset)
  L = Lset(iL);  up=[L,0,0; L+1,0,0];
  for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
  Rr = nan(ndraw,numel(bdb),nS);  PRr = Rr;  PBr = Rr;  Rlin = nan(ndraw,1);  Rbel = nan(ndraw,1);
  for d = 1:ndraw
    rand('state',1000*L+d); randn('state',1000*L+d);
    [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
    crn = make_crn(K,N,S);
    ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax, ...
                  'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
    tic;  gu = ga_design(ctx0,popsize,gens,tsize,mrate);                   % design is independent of beta
    [~,Rbel(d)] = eval_ind(gu,ctx0,false);                                 % rate it believes (linear model)
    printf('L=%d draw=%d GA-unaware designed (%.0fs), linear-model rate %.2f\n',L,d,toc,Rbel(d)); fflush(stdout);
    for ib = 1:numel(bdb)
      A = sqrt(10^(bdb(ib)/10)*Pris/N);
      ctx = ctx0;  ctx.A = A;
      cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
      ctx.mode='linreg'; tic; gl = ga_design(ctx,popsize,gens,tsize,mrate); tl=toc;
      ctx.mode='aware';  tic; ga = ga_design(ctx,popsize,gens,tsize,mrate); ta=toc;
      d0 = cf.aware;  tic; gg = grad_design(ctx,d0,giters); tg=toc;
      D = cell(1,nS);
      D{1} = cf.unaware;  D{2} = cf.aware;
      D{3} = to_design(gu);  D{5} = to_design(ga);  D{6} = to_design(gg);
      % linear-region design: effective gain after the cap
      Pu = sum(abs(G*gl.w).^2,2)+sv2;  D{4} = to_design(gl);  D{4}.p = min(D{4}.p, A*sqrt(rho_max./Pu));
      % unaware designs are re-projected onto the TRUE budget (it can only shrink)
      for s = [1 3]
        ind.w=D{s}.w; ind.P=diag(D{s}.p); ind.Theta=diag(exp(1j*D{s}.theta)); c2=ctx; c2.mode='aware';
        ind = project_ris(ind,c2); D{s}.p = real(diag(ind.P));
      end
      for s = 1:nS
        sinr = sim_sinr(h,G,f,D{s}.theta,D{s}.p,D{s}.w,sv2,s2,A,Tsim,555+d);
        [~,~,PBs_,PRs_] = compute_metrics_nl(h,G,f,D{s}.theta,D{s}.p,D{s}.w,sv2,s2,A,crn);
        Rr(d,ib,s) = sum(log2(1+sinr));  PRr(d,ib,s) = PRs_/Pris;  PBr(d,ib,s) = PBs_/Pbs;
      end
      printf('L=%d draw=%d beta=%+d dB | ',L,d,bdb(ib));
      for s=1:nS, printf('%s %.2f  ',schemes{s},Rr(d,ib,s)); end
      printf('| budget%% max %.0f | t(GA-lin %.0fs GA-aw %.0fs grad %.0fs)\n',100*max(squeeze(PRr(d,ib,:))),tl,ta,tg); fflush(stdout);
    end
    save(sprintf('res_%s_L%d.mat',tag,L),'Rr','PRr','PBr','Rbel','bdb','schemes','L','d');
  end
end
