function sinr = sim_sinr(h, G, f, theta, p, w, sv2, s2, A, T, seed)
% SIM_SINR  Symbol-level "ground truth": Gaussian symbols are passed through the actual Rapp
% amplifier element by element and the end-to-end SINR of every user is estimated from the samples.
%   sinr_k = |a_k|^2 / E|r_k - a_k s_k|^2,  a_k = E[r_k s_k*]
if nargin < 11, seed = 12345; end
if isinf(A), A = 1e30; end
K = size(h,1); N = numel(p);  ph = exp(1j*theta(:));  p = p(:);
st = rand('state'); sn = randn('state'); rand('state',seed); randn('state',seed);
chunk = 10000; nch = ceil(T/chunk);
Sy = zeros(K,0); R = zeros(K,0);
Sall = complex(zeros(K,nch*chunk)); Rall = complex(zeros(K,nch*chunk));
for c = 1:nch
    cn = @(a,b) (randn(a,b)+1j*randn(a,b))/sqrt(2);
    s = cn(K,chunk);  x = w*s;
    u = G*x + sqrt(sv2)*cn(N,chunk);
    y = (p.*u) ./ (1 + (p.*abs(u)/A).^4).^0.25 .* ph;
    r = h*x + f*y + sqrt(s2)*cn(K,chunk);
    idx = (c-1)*chunk+(1:chunk);  Sall(:,idx) = s;  Rall(:,idx) = r;
end
rand('state',st); randn('state',sn);
sinr = zeros(K,1);
for k = 1:K
    a = mean(Rall(k,:).*conj(Sall(k,:)));
    res = Rall(k,:) - a*Sall(k,:);
    sinr(k) = abs(a)^2 / mean(abs(res).^2);
end
end
