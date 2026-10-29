function out = closed_form_design(h, G, f, sv2, s2, Pbs, Pris, A, pmax, crn)
% CLOSED_FORM_DESIGN  Structured reference designs (no iterative optimiser):
%   phases aligned to the cascaded path of user 1, ZF precoder on the equivalent channel,
%   one uniform small-signal gain p.
%     out.unaware : p set by the LINEAR power budget (what a compression-blind design does)
%     out.aware   : p chosen by line search on the compression-aware metric under the true budget
% Each field holds theta, p (N x 1), w and the model-M2 rate / powers.
[K, M] = size(h);  N = size(G,1);
W0 = h';  W0 = W0/norm(W0,'fro')*sqrt(Pbs);                % MRT start, as in Figure12_saturation
th = angle(conj(f(1,:).')) - angle(G*W0(:,1));
zf = @(psi) zf_precoder(h, G, f, psi, Pbs);
% --- compression-unaware: largest p meeting the linear RIS budget
lo = 1; hi = pmax;
pl = @(p) pris_lin(h,G,f,th,p,zf,sv2,s2);
if pl(hi) <= Pris, plin = hi; else
    for it = 1:40
        mid = sqrt(lo*hi);
        if pl(mid) <= Pris, lo = mid; else hi = mid; end
    end
    plin = lo;
end
W = zf(plin*exp(1j*th));
out.unaware.theta = th; out.unaware.p = plin*ones(N,1); out.unaware.w = W;
% --- compression-aware: grid over p, fixed-point precoder on the effective reflection
best = -inf;  grid = logspace(0, log10(pmax), 60);
for p = grid
    W = zf(p*exp(1j*th));
    for it = 1:2
        Pu = sum(abs(G*W).^2, 2) + sv2;
        kap = rapp_bussgang(Pu, p*ones(N,1), A);
        W = zf(kap.*exp(1j*th));
    end
    [~, R, ~, PR] = compute_metrics_nl(h,G,f,th,p*ones(N,1),W,sv2,s2,A,crn);
    if PR <= Pris && R > best
        best = R;  bp = p;  bW = W;
    end
end
if isinf(best), bp = 1; bW = zf(exp(1j*th)); end
out.aware.theta = th; out.aware.p = bp*ones(N,1); out.aware.w = bW;
end

function W = zf_precoder(h, G, f, psi, Pbs)
Heq = h + (f .* psi.') * G;
W = Heq' / (Heq*Heq');
W = W * sqrt(Pbs)/norm(W,'fro');
end

function PR = pris_lin(h,G,f,th,p,zf,sv2,s2)
N = size(G,1);  W = zf(p*exp(1j*th));
[~,~,~,PR] = compute_metrics_nl(h,G,f,th,p*ones(N,1),W,sv2,s2,Inf,[]);
end
