function [best, hist] = ga_design(ctx, popsize, gens, tsize, mrate)
% GA_DESIGN  Joint BS-precoder / RIS-gain / RIS-phase genetic algorithm (adaptive penalty).
% Same chromosome, selection, circular-aware crossover and mutation steps as the IJECE
% package; only the objective model (ctx.mode) differs.  The returned best chromosome is
% projected onto BOTH budgets, so reported rates are for feasible designs.
M=ctx.M; N=ctx.N; K=ctx.K;
pop = cell(popsize,1);
for i = 1:popsize
    ind.Theta = diag(exp(1j*2*pi*rand(N,1)));
    ind.P = diag(1 + (ctx.pmax-1)*rand(N,1));
    w = (randn(M,K)+1j*randn(M,K))/sqrt(2);  ind.w = w/norm(w,'fro')*sqrt(ctx.Pbs);
    pop{i} = project_ris(ind, ctx);
end
lam = 10;  best_fit = -inf;  best = pop{1};  hist = zeros(gens,1);
for g = 1:gens
    fit = zeros(popsize,1);  feas = 0;
    for i = 1:popsize
        [R,PB,PR] = eval_ind(pop{i}, ctx, false);
        vb = max(0, PB/ctx.Pbs - 1);  vr = max(0, PR/ctx.Pris - 1);
        feas = feas + (vb <= 1e-9 && vr <= 1e-9);
        fit(i) = R;  viol(i) = vb + vr;
    end
    if feas/popsize <= 0.3, lam = lam*1.5; else, lam = lam*0.9; end
    fit = fit - lam*viol(:);
    [fb, ib] = max(fit);
    if fb > best_fit, best_fit = fb; best = pop{ib}; end
    hist(g) = best_fit;
    sel = tournament_selection(fit, tsize);
    off = generate_offspring(pop(sel), popsize, @ris_crossover);
    for i = 1:popsize, off{i} = mutate_generic(off{i}, mrate, ctx); end
    [~, iw] = min(fit(sel));  off{iw} = best;
    pop = off;
end
pw = sum(abs(best.w(:)).^2);
if pw > ctx.Pbs, best.w = best.w*sqrt(ctx.Pbs/pw); end
best = project_ris(best, ctx);
end
