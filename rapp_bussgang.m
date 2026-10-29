function [kappa, sd2, Pout] = rapp_bussgang(Pu, p, A)
% RAPP_BUSSGANG  Bussgang gain, distortion power and mean output power of a Rapp
% (smoothness s = 2) soft limiter  g(u) = p*u ./ (1 + (p*|u|/A).^4).^(1/4)
% for a circular complex Gaussian input u ~ CN(0, Pu).  Vectorised over elements.
%   kappa = E[u* g(u)]/Pu,   Pout = E|g(u)|^2,   sd2 = Pout - kappa^2*Pu.
% Expectations over t = |u|^2/Pu ~ Exp(1) use 96-point Gauss-Laguerre quadrature.
persistent T W
if isempty(T)
    n = 96; i = (1:n-1)';
    J = diag(2*(0:n-1)+1) + diag(i,1) + diag(i,-1);     % Laguerre Jacobi matrix
    [V, D] = eig(J);
    [T, o] = sort(diag(D));  w0 = V(1,:)'.^2;  W = w0(o);   % Golub-Welsch weights (mu0 = 1)
end
Pu = Pu(:);  p = p(:);
rho = (p.^2) .* Pu ./ (A.^2);
r   = rho * T.';
kI  = ((1 + r.^2).^(-0.25) .* T.') * W;
oI  = ((1 + r.^2).^(-0.5 ) .* T.') * W;
kappa = p .* kI;
Pout  = (p.^2) .* Pu .* oI;
sd2   = max(Pout - kappa.^2 .* Pu, 0);
end
