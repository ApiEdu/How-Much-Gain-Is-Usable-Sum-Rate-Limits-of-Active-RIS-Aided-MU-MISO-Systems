% VALIDATE_FULL  Model vs symbol-level simulation over many designs.
% env: LV (distance), NCF (CF draws), NGR (gradient draws), GITERS
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; rho_max=0.1; S=160; Tsim=100000;
L=str2num(getenv('LV')); ncf=str2num(getenv('NCF')); ngr=str2num(getenv('NGR')); giters=str2num(getenv('GITERS'));
bdb=[-5 0 5]; up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
for d = 1:max(ncf,ngr)
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  ctx0 = struct('h',h,'G',G,'f',f,'sv2',sv2,'s2',s2,'Pbs',Pbs,'Pris',Pris,'pmax',pmax,'rho_max',rho_max,'M',M,'N',N,'K',K,'crn',crn,'A',Inf,'mode','aware');
  for ib = 1:numel(bdb)
    A = sqrt(10^(bdb(ib)/10)*Pris/N); ctx = ctx0; ctx.A = A;
    cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
    grp = {}; ds = {};
    if d <= ncf, grp{end+1}='CF-aware'; ds{end+1}=cf.aware; grp{end+1}='CF-unaware'; ds{end+1}=cf.unaware; end
    if d <= ngr
      g = grad_design(ctx, cf.aware, giters); grp{end+1}='Grad-aware'; ds{end+1}=to_design(g);
    end
    for q = 1:numel(ds)
      D = ds{q};
      sM = compute_metrics_nl(h,G,f,D.theta,D.p,D.w,sv2,s2,A,crn);
      sI = sinr_indep(h,G,f,D.theta,D.p,D.w,sv2,s2,A);
      sS = sim_sinr(h,G,f,D.theta,D.p,D.w,sv2,s2,A,Tsim,555+d);
      for k=1:K
        printf('VAL %s %d %d %d %d %.4f %.4f %.4f\n',grp{q},L,d,bdb(ib),k,10*log10(sM(k)),10*log10(sI(k)),10*log10(sS(k)));
      end
      fflush(stdout);
    end
  end
end
