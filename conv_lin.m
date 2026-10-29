% CONV_LIN  Same as conv_grad.m but for the linear-region-cap design (equal iteration budget comparison).
% Logs model objective and symbol-level sum-rate at iteration checkpoints. env: LV, GITERS (e.g. 150)
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
  for ib = 1:numel(bdb)
    A = sqrt(10^(bdb(ib)/10)*Pris/N);
    ctx = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',A,'mode','linreg');
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    d0 = cf.aware;  Pu0 = sum(abs(G*d0.w).^2,2)+sv2;  d0.p = max(min(d0.p, A*sqrt(rho_max./Pu0)),1+1e-3);
    tic; [~,hist] = grad_design_h(ctx, d0, giters); tt=toc;
    its = [hist.it];
    for c = cps
      if c > giters, continue; end
      q = find(its==c,1); x = hist(q).x;
      w = reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K); p = 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))); th = x(2*M*K+N+(1:N));
      ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
      pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
      ind = project_ris(ind, ctx);  D = to_design(ind);  Pu = sum(abs(G*D.w).^2,2)+sv2;  D.p = min(D.p, A*sqrt(rho_max./Pu));
      s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
      printf('CONVLIN %d %d %d it %d  J(model, best so far) %.3f  simrate %.3f  (total %.0fs)\n',L,d,bdb(ib),c,hist(q).J,sum(log2(1+s)),tt); fflush(stdout);
    end
  end
end
