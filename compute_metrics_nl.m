function [sinr, R, P_BS, P_RIS] = compute_metrics_nl(h, G, f, theta, p, w, sv2, s2, A, crn)
% COMPUTE_METRICS_NL  Sum-rate of the active-RIS MU-MISO link with a Rapp-type amplifier.
%   h : K x M   f : K x N   G : N x M   theta,p : N x 1   w : M x K   A : saturation amplitude
%   A = Inf reproduces the linear model of compute_metrics2.m exactly.
%
%   Model M2.  Each amplifier is replaced by its Bussgang decomposition
%       y_n = kappa_n*u_n + d_n,   d_n uncorrelated with u.
%   kappa_n and the output power come from quadrature (rapp_bussgang).  The distortion
%   d is NOT assumed independent across elements: its variance at each user is estimated
%   from the common random samples in crn, so the coherent focusing of the distortion
%   by the RIS phases is retained.  (An independent-distortion model was found to be wrong
%   by up to ~4 dB; see validate_m2.m.)
K = size(h,1);  N = numel(p);  theta = theta(:);  p = p(:);
Gw = G*w;
Pu = sum(abs(Gw).^2, 2) + sv2;
ph = exp(1j*theta);
if isinf(A)
    kap = p;  Pout = p.^2 .* Pu;  dist = zeros(K,1);
else
    [kap, ~, Pout] = rapp_bussgang(Pu, p, A);
    if nargin >= 10 && ~isempty(crn)
        U = Gw*crn.s + sqrt(sv2)*crn.v;                            % N x S
        Y = (p.*U) ./ (1 + (p.*abs(U)/A).^4).^0.25;
        D = Y - kap.*U;
        dist = mean(abs(f*(ph.*D)).^2, 2);                         % K x 1
    else
        error('compute_metrics_nl: pass crn for finite A');
    end
end
Heq = h + (f .* (kap.*ph).') * G;                                 % K x M
Sm  = abs(Heq*w).^2;
sig = diag(Sm);  intf = sum(Sm,2) - sig;
noise = (abs(f).^2)*(kap.^2*sv2) + dist + s2;
sinr = sig ./ (intf + noise);
R = sum(log2(1+sinr));
P_BS = sum(abs(w(:)).^2);
P_RIS = sum(Pout);
end
