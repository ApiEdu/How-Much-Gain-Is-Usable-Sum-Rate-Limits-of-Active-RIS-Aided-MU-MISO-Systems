% EVAL_FP  Score the linear-model FP/WMMSE design (J98-style) on the compressed hardware by symbol-level simulation.
more off; sv2=1e-10; s2=1e-10; N=256; Pris=0.1; rho_max=0.1; Tsim=40000;
bdb=[-5 0 5];
for L=[200 100]
  R=nan(14,3,2);   % (draw, beta, [FP as designed, FP + linear-region cap])
  for d=1:14
    c=load(sprintf('ch/ch_L%d_d%d.mat',L,d)); q=load(sprintf('ch/fp_L%d_d%d.mat',L,d));
    p=q.p(:); th=q.theta(:); w=q.w; Pu=sum(abs(c.G*w).^2,2)+sv2;
    for ib=1:3
      A=sqrt(10^(bdb(ib)/10)*Pris/N);
      s=sim_sinr(c.h,c.G,c.f,th,p,w,sv2,s2,A,Tsim,555+d); R(d,ib,1)=sum(log2(1+s));
      pc=min(p, A*sqrt(rho_max./Pu));
      s=sim_sinr(c.h,c.G,c.f,th,pc,w,sv2,s2,A,Tsim,555+d); R(d,ib,2)=sum(log2(1+s));
    end
  end
  save(sprintf('res_fp_L%d.mat',L),'R','bdb');
  printf('L=%d: beta -5/0/5 FP-as-designed: %s | FP+cap: %s\n',L,mat2str(squeeze(mean(R(:,:,1),1)),4),mat2str(squeeze(mean(R(:,:,2),1)),4));
end
