function sinr = sinr_indep(h,G,f,theta,p,w,sv2,s2,A)
% SINR_INDEP  Bussgang model that treats the distortion as independent across elements
% (no coherent focusing by the RIS phases), for comparison with compute_metrics_nl.
K = size(h,1); p = p(:); ph = exp(1j*theta(:));
Pu = sum(abs(G*w).^2,2) + sv2;
[kap, sd2, ~] = rapp_bussgang(Pu, p, A);
dist = (abs(f).^2)*sd2;
Heq = h + (f .* (kap.*ph).') * G;
Sm = abs(Heq*w).^2; sig = diag(Sm); intf = sum(Sm,2)-sig;
noise = (abs(f).^2)*(kap.^2*sv2) + dist + s2;
sinr = sig./(intf+noise);
end
