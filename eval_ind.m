function [R, PB, PR, sinr] = eval_ind(ind, ctx, truth)
% EVAL_IND  Objective seen by the GA for each design philosophy (ctx.mode):
%   'unaware' : GA optimises the LINEAR active-RIS model (as in the IJECE paper);
%               truth=true re-scores the same chromosome on the compressed hardware.
%   'linreg'  : linear-region restriction (Yang-style): gain is capped so that the average
%               element stays below rho_max*Asat^2 at its input, then scored on the true model.
%   'aware'   : compression-aware (model M2 throughout).
th = angle(diag(ind.Theta));  p = real(diag(ind.P));  w = ind.w;
if strcmp(ctx.mode,'unaware') && ~truth, Ae = Inf; else, Ae = ctx.A; end
if strcmp(ctx.mode,'linreg')
    Pu = sum(abs(ctx.G*w).^2, 2) + ctx.sv2;
    p = min(p, ctx.A*sqrt(ctx.rho_max./Pu));
end
[sinr, R, PB, PR] = compute_metrics_nl(ctx.h, ctx.G, ctx.f, th, p, w, ctx.sv2, ctx.s2, Ae, ctx.crn);
end
