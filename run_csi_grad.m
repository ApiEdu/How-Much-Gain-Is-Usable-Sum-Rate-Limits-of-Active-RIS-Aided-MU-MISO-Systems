% RUN_CSI_GRAD  Linear-model and capped gradient designs computed from imperfect channel estimates (NMSE -10 dB, independent
% complex Gaussian errors on h, G and f, relative to the mean channel power of each matrix), scored on the TRUE channels.
% 150 iterations, beta = 0 dB, draws 1..4 at distance LV. Output lines: CSIG L d nmse beta design iters rate
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000; nm=-10; b=0; giters=str2num(getenv('GITERS'));
L=str2num(getenv('LV')); up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
cn = @(sz) (randn(sz)+1j*randn(sz))/sqrt(2);
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
A = sqrt(10^(b/10)*Pris/N); e = 10^(nm/10);
for d = 1:4
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  randn('state',9999+1000*L+d);
  hh = h + sqrt(e*mean(abs(h(:)).^2))*cn(size(h));
  Gh = G + sqrt(e*mean(abs(G(:)).^2))*cn(size(G));
  fh = f + sqrt(e*mean(abs(f(:)).^2))*cn(size(f));
  ctx0 = struct('h',hh,'G',Gh,'f',fh,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
  % --- linear-model design from the estimates, scored on the true channels at beta
  cf0 = closed_form_design(hh,Gh,fh,sv2,s2,Pbs,Pris,1e3,pmax,crn);
  [~,hist] = grad_design_h(ctx0, cf0.unaware, giters); its=[hist.it];
  q = find(its==giters,1); [w,p,th] = dec(hist(q).x);
  ctx = ctx0; ctx.A = A; ctx.mode='aware';
  ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
  pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
  ind = project_ris(ind, ctx); D = to_design(ind);
  s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
  printf('CSIG %d %d %d %d unaware %d %.3f\n',L,d,nm,b,giters,sum(log2(1+s))); fflush(stdout);
  % --- capped design from the estimates
  cf = closed_form_design(hh,Gh,fh,sv2,s2,Pbs,Pris,A,pmax,crn);
  ctx = ctx0; ctx.A = A; ctx.mode = 'linreg'; d0 = cf.aware;
  Pu0 = sum(abs(Gh*d0.w).^2,2)+sv2; d0.p = max(min(d0.p, A*sqrt(rho_max./Pu0)),1+1e-3);
  [~,hist] = grad_design_h(ctx, d0, giters); its=[hist.it];
  q = find(its==giters,1); [w,p,th] = dec(hist(q).x);
  ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
  pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
  ind = project_ris(ind, ctx); D = to_design(ind);
  Pu = sum(abs(Gh*D.w).^2,2)+sv2; D.p = min(D.p, A*sqrt(rho_max./Pu));
  s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
  printf('CSIG %d %d %d %d linreg %d %.3f\n',L,d,nm,b,giters,sum(log2(1+s))); fflush(stdout);
end
