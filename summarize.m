% SUMMARIZE  Mean +/- std over channel draws of the symbol-level sum-rate, per L, beta, design.
more off;
for L = [200 100]
  load(sprintf('res_main_L%d.mat',L));       % Rr(draw,beta,scheme), PRr, PBr, Rbel, bdb, schemes
  nd = size(Rr,1);
  printf('\n=== L = %d m  (%d channel draws; symbol-level sum-rate, bit/s/Hz) ===\n', L, nd);
  printf('%-8s', 'beta dB'); printf('%-16s', schemes{:}); printf('\n');
  for ib = 1:numel(bdb)
    printf('%-8d', bdb(ib));
    for s = 1:numel(schemes), x = Rr(:,ib,s); printf('%5.2f+-%-8.2f', mean(x), std(x)); end
    printf('\n');
  end
  printf('GA belief under the linear model (compression ignored): %.2f +- %.2f\n', mean(Rbel), std(Rbel));
  printf('RIS budget use (mean %% of P_RIS,max):\n');
  for ib = 1:numel(bdb)
    printf('%-8d', bdb(ib)); for s = 1:numel(schemes), printf('%-16.0f', 100*mean(PRr(:,ib,s))); end; printf('\n');
  end
  printf('Paired differences (mean over draws):\n');
  for ib = 1:numel(bdb)
    d1 = Rr(:,ib,5)-Rr(:,ib,4); d2 = Rr(:,ib,5)-Rr(:,ib,3); d3 = Rr(:,ib,6)-Rr(:,ib,5); d4 = Rr(:,ib,6)-Rr(:,ib,2); d5 = Rr(:,ib,4)-Rr(:,ib,3);
    printf('beta %+3d dB: GA-aware - GA-linregion %+6.2f (draws>0: %d/%d) | GA-aware - GA-unaware %+6.2f | GA-linregion - GA-unaware %+6.2f | Grad - GA-aware %+6.2f | Grad - CF-aware %+6.2f\n', ...
           bdb(ib), mean(d1), sum(d1>0), nd, mean(d2), mean(d5), mean(d3), mean(d4));
  end
end
