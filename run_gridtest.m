% Is the small deficit of the aware closed-form design at high beta a p-grid effect? 60-point vs 480-point grid, beta=30 dB
more off;
M=4; N=256; K=2; total=10; Pbs=0.99*total; Pris=0.01*total; sv2=1e-10; s2=1e-10;
fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
pmax=100; S=160; Tsim=40000; b=30;
L=100; up=[L,0,0; L+1,0,0];
for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
for d=1:3
  rand('state',1000*L+d); randn('state',1000*L+d);
  [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
  crn = make_crn(K,N,S);
  A = sqrt(10^(b/10)*Pris/N);
  c1 = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
  c2 = closed_form_fine(h,G,f,sv2,s2,Pbs,Pris,A,pmax,crn);
  su=sim_sinr(h,G,f,c1.unaware.theta,c1.unaware.p,c1.unaware.w,sv2,s2,A,Tsim,555+d);
  s1=sim_sinr(h,G,f,c1.aware.theta,c1.aware.p,c1.aware.w,sv2,s2,A,Tsim,555+d);
  s2_=sim_sinr(h,G,f,c2.aware.theta,c2.aware.p,c2.aware.w,sv2,s2,A,Tsim,555+d);
  printf('GRID %d %d unaware %.3f aware60 %.3f aware480 %.3f  p: unaware %.3f aware60 %.3f aware480 %.3f\n',L,d,sum(log2(1+su)),sum(log2(1+s1)),sum(log2(1+s2_)),c1.unaware.p(1),c1.aware.p(1),c2.aware.p(1)); fflush(stdout);
end
