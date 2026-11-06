% RUN_LONG  Long-budget check (K=2, N=256): capped (linreg) and compression-aware (aware) gradient designs, each started from
% its own closed-form design (as in run_main150.m), one channel draw DRAW at distance LV, saturation ratios BETAS (dB),
% scored at the checkpoints CPS (multiples of 10, last = GITERS). Output lines start with LONG:
%   LONG <L> <draw> <beta_dB> <mode> <iteration> <symbol-level sum-rate> <design-time objective J>
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000;
L=str2num(getenv('LV')); d=str2num(getenv('DRAW')); giters=str2num(getenv('GITERS'));
bdb=str2num(getenv('BETAS')); cps=str2num(getenv('CPS')); up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
rand('state',1000*L+d); randn('state',1000*L+d);
[h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
crn = make_crn(K,N,S);
ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
for b = bdb
  A = sqrt(10^(b/10)*Pris/N);  base = ctx0; base.A = A;
  cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
  for mode = {'aware','linreg'}
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
      printf('LONG %d %d %d %s %d %.3f %.3f\n',L,d,b,mode{1},c,sum(log2(1+s)),hist(q).J); fflush(stdout);
    end
  end
end
