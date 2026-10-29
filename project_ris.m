function ind = project_ris(ind, ctx)
% PROJECT_RIS  Shrink the amplification vector (exact bisection on a common scale) until the
% RIS power budget holds under the model of ctx.mode.  No-op when already feasible.
p0 = real(diag(ind.P));
ind.P = diag(min(max(p0,1),ctx.pmax));  p0 = real(diag(ind.P));
[~,~,PR] = eval_ind(ind, ctx, false);
if PR <= ctx.Pris, return; end
lo = 0; hi = 1;
for it = 1:28
    s = (lo+hi)/2;
    ind.P = diag(max(s*p0,1));
    [~,~,PR] = eval_ind(ind, ctx, false);
    if PR <= ctx.Pris, lo = s; else, hi = s; end
end
ind.P = diag(max(lo*p0,1));
end
