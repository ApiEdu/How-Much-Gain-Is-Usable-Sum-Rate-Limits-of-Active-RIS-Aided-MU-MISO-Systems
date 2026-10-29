M=4; N=256; K=2; fc=5; nAI=5; nDI=2; BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
for L=[200 100]
  up=[L,0,0; L+1,0,0];
  for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
  for d=1:14
    rand('state',1000*L+d); randn('state',1000*L+d);
    [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up); h=h'; f=f';
    save('-v6',sprintf('ch/ch_L%d_d%d.mat',L,d),'h','f','G');
  end
end
