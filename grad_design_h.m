function [best, hist] = grad_design_h(ctx, d0, iters)
% GRAD_DESIGN  Projected gradient ascent (Adam, finite-difference gradients) on the
% compression-aware objective, warm-started from design d0 (struct with theta,p,w).
% Used as the "smooth-optimiser" benchmark that a reviewer will ask for when GA is proposed.
M=ctx.M; N=ctx.N; K=ctx.K; pmax=ctx.pmax;
ctxg = ctx;  ctxg.crn.s = ctx.crn.s(:,1:96);  ctxg.crn.v = ctx.crn.v(:,1:96);   % cheaper gradient
dec = @(x) deal(reshape(x(1:M*K)+1j*x(M*K+1:2*M*K),M,K), 1+(pmax-1)./(1+exp(-x(2*M*K+(1:N)))), x(2*M*K+N+(1:N)));
p0 = min(max(d0.p,1+1e-3),pmax-1e-3);
x = [real(d0.w(:)); imag(d0.w(:)); log((p0-1)./(pmax-p0)); d0.theta(:)];
nv = numel(x);
lr = [0.03*ones(2*M*K,1); 0.2*ones(N,1); 0.05*ones(N,1)];
st = [1e-4*ones(2*M*K,1); 1e-3*ones(N,1); 1e-4*ones(N,1)];
m = zeros(nv,1); v = zeros(nv,1);
J = @(x,c) objective(x,c,dec);
bestJ = -inf; bestx = x; hist = struct('it',{},'J',{},'x',{});
for it = 1:iters
    J0 = J(x,ctxg);  gr = zeros(nv,1);
    for i = 1:nv
        xi = x; xi(i) = xi(i) + st(i);
        gr(i) = (J(xi,ctxg) - J0)/st(i);
    end
    m = 0.9*m + 0.1*gr;  v = 0.999*v + 0.001*gr.^2;
    x = x + lr.*(m/(1-0.9^it))./(sqrt(v/(1-0.999^it)) + 1e-8);
    if mod(it,10) == 0 || it == iters
        Jx = J(x,ctx);
        if Jx > bestJ, bestJ = Jx; bestx = x; end
        hist(end+1) = struct('it',it,'J',bestJ,'x',bestx);
    end
end
[w,p,th] = dec(bestx);
ind.w = w; ind.P = diag(p); ind.Theta = diag(exp(1j*th(:)));
pw = sum(abs(ind.w(:)).^2);  if pw > ctx.Pbs, ind.w = ind.w*sqrt(ctx.Pbs/pw); end
best = project_ris(ind, ctx);
end

function j = objective(x, ctx, dec)
[w,p,th] = dec(x);
if strcmp(ctx.mode,'linreg'), Pu = sum(abs(ctx.G*w).^2,2)+ctx.sv2; p = min(p, ctx.A*sqrt(ctx.rho_max./Pu)); end
[~,R,PB,PR] = compute_metrics_nl(ctx.h,ctx.G,ctx.f,th(:),p,w,ctx.sv2,ctx.s2,ctx.A,ctx.crn);
j = R - 30*max(0,PR/ctx.Pris-1) - 30*max(0,PB/ctx.Pbs-1);
end
