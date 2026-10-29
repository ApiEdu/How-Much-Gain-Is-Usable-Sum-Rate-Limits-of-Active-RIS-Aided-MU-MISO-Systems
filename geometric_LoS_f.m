function H_LoS = geometric_LoS_f(N, RIS_pos, user_pos, lambda)
% Geometric LoS component for RIS-User channel
N_sqrt = round(sqrt(N));
d_RIS = lambda/2;
H_LoS = zeros(N, 1);

for n = 1:N
    row = mod(n-1, N_sqrt);
    col = floor((n-1)/N_sqrt);
    RIS_ant_pos = RIS_pos + [row*d_RIS, col*d_RIS, 0];
    
    distance = norm(user_pos - RIS_ant_pos);
    phase = 2*pi*distance/lambda;
    H_LoS(n) = exp(-1j*phase);
end

H_LoS = H_LoS / sqrt(N);
end