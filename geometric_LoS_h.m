function H_LoS = geometric_LoS_h(M, BS_pos, user_pos, lambda)
% Geometric LoS component for BS-User channel
d_BS = lambda/2;
H_LoS = zeros(M, 1);

for m = 1:M
    BS_ant_pos = BS_pos + [0, (m-1)*d_BS, 0];
    distance = norm(user_pos - BS_ant_pos);
    phase = 2*pi*distance/lambda;
    H_LoS(m) = exp(-1j*phase);
end

H_LoS = H_LoS / sqrt(M);
end