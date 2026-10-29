% RUN_MAIN150  Main comparison at a converged-ish budget: linear-model (unaware), capped (linreg) and compression-aware
% gradient designs, each GITERS iterations, scored at checkpoints 50 and GITERS. Draws 1..4 at distance LV.
% Unaware design is run once per draw (independent of beta) and scored at beta -5/0/5; aware and linreg run only for the
% (draw:beta) pairs in PAIRS, e.g. "1:0;2:0;3:-5". Output lines start with HBG.
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000;
L=str2num(getenv('LV')); giters=str2num(getenv('GITERS')); up=[L,0,0; L+1,0,0];
pairs=strsplit(getenv('PAIRS'),';'); pd=zeros(numel(pairs),1); pb=pd;
for i=1:numel(pairs), q=strsplit(pairs{i},':'); pd(i)=str2num(q{1}); pb(i)=str2num(q{2}); end
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
cps=[50 giters]; allb=[10 20];
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
for d = 1:4
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
  % --- unaware (once per draw)
  cf0 = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,1e3,pmax,crn);
  [~,hist] = grad_design_h(ctx0, cf0.unaware, giters); its=[hist.it];
  for c = cps
    q = find(its==c,1); [w,p,th] = dec(hist(q).x);
    for ib = 1:numel(allb)
      A = sqrt(10^(allb(ib)/10)*Pris/N); ctx = ctx0; ctx.A = A; ctx.mode='aware';
      ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
      pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
      ind = project_ris(ind, ctx); D = to_design(ind);
      s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
      printf('HBG %d %d %d unaware %d %.3f\n',L,d,allb(ib),c,sum(log2(1+s))); fflush(stdout);
    end
  end
  % --- aware and capped for the requested pairs
  for i = find(pd==d)'
    b = pb(i); A = sqrt(10^(b/10)*Pris/N);
    base = ctx0; base.A = A;
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    for mode = {'linreg'}
      ctx = base; ctx.mode = mode{1}; d0 = cf.aware;
      if strcmp(mode{1},'linreg'), Pu0 = sum(abs(G*d0.w).^2,2)+sv2; d0.p = max(min(d0.p, A*sqrt(rho_max./Pu0)),1+1e-3); end
      [~,hist] = grad_design_h(ctx, d0, giters); its=[hist.it];
      for c = cps
        q = find(its==c,1); [w,p,th] = dec(hist(q).x);
        ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
        pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
        ind = project_ris(ind, ctx); D = to_design(ind);
        if strcmp(mode{1},'linreg'), Pu = sum(abs(G*D.w).^2,2)+sv2; D.p = min(D.p, A*sqrt(rho_max./Pu)); end
        s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
        printf('HBG %d %d %d %s %d %.3f\n',L,d,b,mode{1},c,sum(log2(1+s))); fflush(stdout);
      end
    end
  end
end
