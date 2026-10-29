function s = to_design(ind)
% Chromosome -> plain design struct (theta, p, w)
s.theta = angle(diag(ind.Theta));  s.p = real(diag(ind.P));  s.w = ind.w;
end
