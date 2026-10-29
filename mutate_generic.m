function ind = mutate_generic(ind, rate, ctx)
% Mutation shared by all GA variants (same step sizes as mutate_individual_fixed.m);
% only the power model used by the projection differs with ctx.mode.
M = ctx.M; N = ctx.N; K = ctx.K;
if rand < rate
    ind.w = ind.w + 0.1*(randn(M,K) + 1j*randn(M,K));
    pw = sum(abs(ind.w(:)).^2);
    if pw > ctx.Pbs, ind.w = ind.w*sqrt(ctx.Pbs/pw); end
    ind = project_ris(ind, ctx);
end
if rand < rate
    p = real(diag(ind.P)).*exp(0.2*randn(N,1));
    ind.P = diag(min(max(p,1),ctx.pmax));
    ind = project_ris(ind, ctx);
end
if rand < rate
    ind.Theta = diag(exp(1j*mod(angle(diag(ind.Theta)) + 0.1*randn(N,1), 2*pi)));
end
end
