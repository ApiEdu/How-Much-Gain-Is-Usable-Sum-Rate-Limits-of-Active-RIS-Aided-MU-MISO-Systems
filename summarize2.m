more off;
for L=[200 100]
  m=load(sprintf('res_main_L%d.mat',L)); e=load(sprintf('res_extra_L%d.mat',L));
  Ga=m.Rr(:,:,6); nd=size(Ga,1);
  printf('\n=== L=%d (%d draws) mean sum-rate, bit/s/Hz ===\n',L,nd);
  printf('beta | Grad-aware | Grad-linreg | Grad-unaware | aware-linreg (wins) | aware-unaware\n');
  for ib=1:numel(m.bdb)
    a=Ga(:,ib); l=e.Rlr(:,ib); u=e.Run(:,ib);
    printf('%+3d  | %6.2f | %6.2f | %6.2f | %+6.2f (%d/%d) | %+6.2f\n',m.bdb(ib),mean(a),mean(l),mean(u),mean(a-l),sum(a>l),nd,mean(a-u));
  end
end
