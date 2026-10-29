% EVAL_FPCAP  Score the in-loop-cap WMMSE designs (fp_cap.py) by symbol-level simulation.
more off; sv2=1e-10; s2=1e-10; Tsim=40000; ITERS=200; bdb=[-5 0 5];
for L=[100 200]
  for d=1:4
    c=load(sprintf('ch/ch_L%d_d%d.mat',L,d));
    for ib=1:3
      fn=sprintf('ch/fpcap%d_L%d_d%d_b%d.mat',ITERS,L,d,bdb(ib)+5);
      if ~exist(fn,'file'), continue; end
      q=load(fn); A=sqrt(10^(bdb(ib)/10)*0.1/256);
      s=sim_sinr(c.h,c.G,c.f,q.theta(:),q.p(:),q.w,sv2,s2,A,Tsim,555+d);
      printf('FPCAP %d %d %d simrate %.3f linrate_start %.3f linrate_end %.3f\n',L,d,bdb(ib),sum(log2(1+s)),q.hist(1),q.hist(end)); fflush(stdout);
    end
  end
end
