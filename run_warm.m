% RUN_WARM  "Release-the-cap" test.  From ONE shared starting point (the capped design after CITERS
% iterations, identical to the 'linreg' design of run_main150x.m) two continuations of XITERS
% iterations are run with the same optimiser, the same fresh Adam state and the same budget:
%   capcont   : the linear-region cap is kept
%   awarewarm : the cap is released (compression-aware objective, same RIS power budget)
% The start is a candidate in both, so awarewarm cannot have a lower design-time objective than the start.
% Output lines start with WARM; PST lines give gain statistics of the final design.
1;
function printdj(tag,it,D,h,G,f,sv2,s2,A,crn,crnF,Pbs,Pris,L,d,b)
% DJ: design-time objective on the S=160 design samples, model rate on 2000 fresh samples, RIS power / budget
[~,R1,PB,PR] = compute_metrics_nl(h,G,f,D.theta,D.p,D.w,sv2,s2,A,crn);
J1 = R1 - 30*max(0,PR/Pris-1) - 30*max(0,PB/Pbs-1);
[~,R2] = compute_metrics_nl(h,G,f,D.theta,D.p,D.w,sv2,s2,A,crnF);
printf('DJ %d %d %d %s %d %.3f %.3f %.3f\n',L,d,b,tag,it,J1,R2,PR/Pris); fflush(stdout);
end

more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000;
L=str2num(getenv('LV')); citers=str2num(getenv('CITERS')); xiters=str2num(getenv('XITERS')); up=[L,0,0; L+1,0,0];
DR=str2num(getenv('DRAWS')); BETAS=str2num(getenv('BETAS'));
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
for d = DR
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  st_=rand('state'); sn_=randn('state'); rand('state',777+d); randn('state',777+d); crnF = make_crn(K,N,2000); rand('state',st_); randn('state',sn_);
  ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
  for b = BETAS
    A = sqrt(10^(b/10)*Pris/N);
    base = ctx0; base.A = A;
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    ctxc = base; ctxc.mode = 'linreg';
    ctxa = base; ctxa.mode = 'aware';
    d0 = cf.aware; Pu0 = sum(abs(G*d0.w).^2,2)+sv2; d0.p = max(min(d0.p, A*sqrt(rho_max./Pu0)),1+1e-3);
    % --- stage 1: capped design, as in run_main150x.m (reproduces its 'linreg' rows)
    [~,hist] = grad_design_h(ctxc, d0, citers); its=[hist.it]; q = find(its==citers,1);
    [w,p,th] = dec(hist(q).x);
    ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
    pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
    ind = project_ris(ind, ctxc); D = to_design(ind);
    Pu = sum(abs(G*D.w).^2,2)+sv2; D.p = min(D.p, A*sqrt(rho_max./Pu));
    s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
    printf('WARM %d %d %d cap %d %.3f\n',L,d,b,citers,sum(log2(1+s)));
    printf('PST %d %d %d cap %d %.2f %.2f %.3f\n',L,d,b,citers,mean(D.p),max(D.p),mean(D.p>0.99*pmax)); fflush(stdout);
    start = D;
    printdj('cap',citers,D,h,G,f,sv2,s2,A,crn,crnF,Pbs,Pris,L,d,b);
    % --- stage 2: two continuations from the same start
    for mode = {'capcont','awarewarm'}
      if strcmp(mode{1},'capcont'), ctx = ctxc; else, ctx = ctxa; end
      [~,hist] = grad_design_w(ctx, start, xiters); its=[hist.it];
      for c = unique([round(xiters/2) xiters])
        q = find(its==c,1); [w,p,th] = dec(hist(q).x);
        ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
        pw = sum(abs(ind.w(:)).^2);  if pw > Pbs, ind.w = ind.w*sqrt(Pbs/pw); end
        ind = project_ris(ind, ctx); D = to_design(ind);
        if strcmp(mode{1},'capcont'), Pu = sum(abs(G*D.w).^2,2)+sv2; D.p = min(D.p, A*sqrt(rho_max./Pu)); end
        s = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
        printf('WARM %d %d %d %s %d %.3f\n',L,d,b,mode{1},c,sum(log2(1+s)));
        printf('PST %d %d %d %s %d %.2f %.2f %.3f\n',L,d,b,mode{1},c,mean(D.p),max(D.p),mean(D.p>0.99*pmax)); fflush(stdout);
        printdj(mode{1},c,D,h,G,f,sv2,s2,A,crn,crnF,Pbs,Pris,L,d,b);
      end
    end
  end
end
