% CONV_UN  Linear-model (unaware) gradient design run for GITERS iterations, re-projected on the true budget and scored.
% The design does not depend on beta, so one run per draw serves all beta. env: LV, GITERS
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000;
L=str2num(getenv('LV')); giters=str2num(getenv('GITERS')); bdb=[-5 5]; up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
cps=[50 100 150];
for d = 1:2
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
  cf0 = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,1e3,pmax,crn);
  tic; [~,hist] = grad_design_h(ctx0, cf0.unaware, giters); tt=toc; its=[hist.it];
  for c = cps
    q = find(its==c,1); x = hist(q).x;
    w = reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K); p = 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))); th = x(2*M*K+N+(1:N));
    for ib = 1:numel(bdb)
      A = sqrt(10^(bdb(ib)/10)*Pris/N); ctx = ctx0; ctx.A = A; ctx.mode='aware';
      ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
      pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
      ind = project_ris(ind, ctx); D = to_design(ind);
      s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
      printf('CONVUN %d %d %d it %d simrate %.3f (%.0fs)\n',L,d,bdb(ib),c,sum(log2(1+s)),tt); fflush(stdout);
    end
  end
end
