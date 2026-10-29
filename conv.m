for L=[200 100]
 for t={'main','extra','fp','more','more2','more3','fp6k'}
  fn=sprintf('res_%s_L%d.mat',t{1},L);
  if exist(fn,'file'), S=load(fn); save('-v6',['v6/' fn],'-struct','S'); end
 end
end
