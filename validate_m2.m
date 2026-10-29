% VALIDATE_M2  Checks (1) the Bussgang quadrature against Monte Carlo, (2) the linear limit
% against the original compute_metrics2.m, (3) model M2 against symbol-level simulation.
more off; rng(3);
M=4; N=256; K=2; Pbs=9.9; Pris=0.1; sv2=1e-10; s2=1e-10; fc=5; nAI=5; nDI=2;
BS=[0,-40,0]; RIS=[200,30,0]; dBR=norm(BS-RIS);
printf('--- (1) Bussgang quadrature vs Monte Carlo ---\n');
for c = [1 1 1; 1 3 1; 0.2 5 2]'
    Pu=c(1); p=c(2); A=c(3);
    u = sqrt(Pu/2)*(randn(2e6,1)+1j*randn(2e6,1));  g = p*u./(1+(p*abs(u)/A).^4).^0.25;
    k = real(mean(conj(u).*g))/Pu;  P = mean(abs(g).^2);
    [k2,sd2,P2] = rapp_bussgang(Pu,p,A);
    printf('Pu=%g p=%g A=%g : kappa MC %.5f quad %.5f | Pout MC %.5f quad %.5f\n',Pu,p,A,k,k2,P,P2);
end
for L = [200 100]
    up = [L,0,0; L+1,0,0];
    for k=1:K, dBU(k,1)=norm(BS-up(k,:)); dRU(k,1)=norm(RIS-up(k,:)); end
    [h,f,G] = Channel_generate_geometric(K,N,M,nAI,nDI,dBR,dBU,dRU,fc,BS,RIS,up);  h=h'; f=f';
    crn = make_crn(K,N,256);
    printf('\n=== L = %d m ===\n', L);
    % (2) linear limit
    W = (randn(M,K)+1j*randn(M,K)); W = W/norm(W,'fro')*sqrt(Pbs);  th = 2*pi*rand(N,1); p = 40*ones(N,1);
    [~,R1] = compute_metrics_nl(h,G,f,th,p,W,sv2,s2,Inf,[]);
    [~,R2] = compute_metrics2(h,G,f,diag(exp(1j*th)),diag(p),W,sv2,s2);
    printf('(2) linear limit: new %.6f vs compute_metrics2 %.6f\n', R1, R2);
    % (3) M2 vs simulation
    for bdB = [-5 0 5]
        A = sqrt(10^(bdB/10)*Pris/N);
        cf = closed_form_design(h,G,f,sv2,s2,Pbs,Pris,A,100,crn);
        d = cf.aware;
        sM = compute_metrics_nl(h,G,f,d.theta,d.p,d.w,sv2,s2,A,crn);
        sS = sim_sinr(h,G,f,d.theta,d.p,d.w,sv2,s2,A,60000,7);
        printf('(3) beta %+d dB  p=%.1f | SINR dB: M2 [%.2f %.2f]  simulated [%.2f %.2f] | rate M2 %.2f sim %.2f\n', ...
               bdB, d.p(1), 10*log10(sM), 10*log10(sS), sum(log2(1+sM)), sum(log2(1+sS)));
    end
end
