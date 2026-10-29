% RUN_CSI_CF  Closed-form designs computed from imperfect channel estimates, scored on the true channels (symbol-level simulation).
% Estimate = true + CN(0, nmse * mean|true|^2) per entry (h, G and f), independent errors. Output: CSI L d nmse_dB beta unaware aware
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; S=160; Tsim=40000; bdb=[-5 0 5]; nm=[-20 -10]; ndraw=14;
L=str2num(getenv('LV')); up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
cn = @(sz) (randn(sz)+1j*randn(sz))/sqrt(2);
for d=1:ndraw
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  randn('state',7777+1000*L+d);                                  % estimation errors: separate stream
  for inm=1:numel(nm)
    e = 10^(nm(inm)/10);
    hh = h + sqrt(e*mean(abs(h(:)).^2))*cn(size(h));
    Gh = G + sqrt(e*mean(abs(G(:)).^2))*cn(size(G));
    fh = f + sqrt(e*mean(abs(f(:)).^2))*cn(size(f));
    for ib=1:numel(bdb)
      A = sqrt(10^(bdb(ib)/10)*Pris/N);
      ctx.h=hh; ctx.G=Gh; ctx.f=fh; ctx.sv2=sv2; ctx.s2=s2; ctx.Pbs=Pbs; ctx.Pris=Pris; ctx.pmax=pmax; ctx.A=A; ctx.crn=crn; ctx.mode='aware'; ctx.M=M; ctx.N=N; ctx.K=K;
      cf = closed_form_design(hh,Gh,fh,sv2,s2,Pbs,Pris,A,pmax,crn);
      ind.w=cf.unaware.w; ind.P=diag(cf.unaware.p); ind.Theta=diag(exp(1j*cf.unaware.theta)); ind=project_ris(ind,ctx); Du=to_design(ind);
      su=sim_sinr(h,G,f,Du.theta,Du.p,Du.w,sv2,s2,A,Tsim,555+d);
      sa=sim_sinr(h,G,f,cf.aware.theta,cf.aware.p,cf.aware.w,sv2,s2,A,Tsim,555+d);
      printf('CSI %d %d %d %d %.3f %.3f\n',L,d,nm(inm),bdb(ib),sum(log2(1+su)),sum(log2(1+sa))); fflush(stdout);
    end
  end
end
