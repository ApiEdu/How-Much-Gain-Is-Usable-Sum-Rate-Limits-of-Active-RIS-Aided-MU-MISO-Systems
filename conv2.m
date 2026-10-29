for L=[100 200], for K=[2 3]
  fn=sprintf('res_nk_L%d_N128_K%d.mat',L,K); S=load(fn); save('-v6',['v6/' fn],'-struct','S');
end, end
