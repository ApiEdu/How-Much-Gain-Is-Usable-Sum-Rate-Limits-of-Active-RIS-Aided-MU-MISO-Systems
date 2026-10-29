% RUN_MORE  Additional channel draws (D0+1..D1): Grad-aware(CF start), Grad-linregion, Grad-aware2 (start from linreg solution), Grad-unaware.
more off;
M=4; N=str2num(getenv('NN')); K=str2num(getenv('KK')); total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=40000; giters=str2num(getenv('GITERS'));
L=str2num(getenv('LSET')); bdb=str2num(getenv('BDB')); D0=str2num(getenv('D0')); D1=str2num(getenv('D1')); tag=getenv('TAG');
up=[L+(0:K-1)',zeros(K,2)];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
nd=D1-D0; R=nan(nd,numel(bdb),4);   % 1 aware(CF) 2 linreg 3 aware2 4 unaware
for dd = 1:nd
  d = D0+dd;
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax, ...
                'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','unaware');
  cf0 = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,1e3,pmax,crn);
  gu = grad_design(ctx0, cf0.unaware, giters);
  for ib = 1:numel(bdb)
    A = sqrt(10^(bdb(ib)/10)*Pris/N);  ctx = ctx0;  ctx.A = A;
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    ctx.mode='aware';  ga = grad_design(ctx, cf.aware, giters);  Da = to_design(ga);
    ctx.mode='linreg';
    d0 = cf.aware;  Pu = sum(abs(G*d0.w).^2,2)+sv2;  d0.p = max(min(d0.p, A*sqrt(rho_max./Pu)),1+1e-3);
    gl = grad_design(ctx, d0, giters);  Dl = to_design(gl);  Pu = sum(abs(G*Dl.w).^2,2)+sv2;  Dl.p = min(Dl.p, A*sqrt(rho_max./Pu));
    ctx.mode='aware';  ga2 = grad_design(ctx, Dl, giters);  Da2 = to_design(ga2);
    ind = project_ris(gu, ctx);  Du = to_design(ind);
    Ds = {Da, Dl, Da2, Du};
    for s=1:4
      sinr = sim_sinr(h,G,f,Ds{s}.theta,Ds{s}.p,Ds{s}.w,sv2,s2,A,Tsim,555+d);  R(dd,ib,s)=sum(log2(1+sinr));
    end
    printf('L=%d draw=%d beta=%+d | aware %.2f linreg %.2f aware2 %.2f unaware %.2f\n',L,d,bdb(ib),R(dd,ib,:)); fflush(stdout);
  end
  save(sprintf('res_%s_L%d_N%d_K%d.mat',tag,L,N,K),'R','bdb','L','D0','D1','dd');
end
