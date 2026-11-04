% RUN_RHO  Sensitivity of the capped design to the linear-region ratio rho_max (env RHOS), same optimiser,
% start rule and budget as the 'linreg' rows of run_main150x.m (rho_max = 0.1 there).
% Env PMAX sets the upper gain bound (100 in the paper).  Output lines start with RHO: L draw beta rho pmax iters sum-rate.   PST lines: gain statistics.
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=str2num(getenv('PMAX')); S=160; Tsim=40000;
L=str2num(getenv('LV')); giters=str2num(getenv('GITERS')); up=[L,0,0; L+1,0,0];
DR=str2num(getenv('DRAWS')); BETAS=str2num(getenv('BETAS')); RHOS=str2num(getenv('RHOS'));
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
for d = DR
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  for b = BETAS
    A = sqrt(10^(b/10)*Pris/N);
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    for rho = RHOS
      ctx = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho,'M',M,'N',N,'K',K,'crn',crn,'A',A,'mode','linreg');
      d0 = cf.aware; Pu0 = sum(abs(G*d0.w).^2,2)+sv2; d0.p = max(min(d0.p, A*sqrt(rho./Pu0)),1+1e-3);
      [~,hist] = grad_design_h(ctx, d0, giters); its=[hist.it]; q = find(its==giters,1);
      [w,p,th] = dec(hist(q).x);
      ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
      pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
      ind = project_ris(ind, ctx); D = to_design(ind);
      Pu = sum(abs(G*D.w).^2,2)+sv2; D.p = min(D.p, A*sqrt(rho./Pu));
      s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
      printf('RHO %d %d %d %.3f %d %d %.3f\n',L,d,b,rho,pmax,giters,sum(log2(1+s)));
      printf('PST %d %d %d rho%.3f %d %.2f %.2f %.3f\n',L,d,b,rho,giters,mean(D.p),max(D.p),mean(D.p>0.99*pmax)); fflush(stdout);
    end
  end
end
