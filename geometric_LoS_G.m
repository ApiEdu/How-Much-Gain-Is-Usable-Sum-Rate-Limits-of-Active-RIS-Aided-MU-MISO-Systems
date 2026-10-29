function H_LoS = geometric_LoS_G(M, N, BS_pos, RIS_pos, lambda)
% Geometric LoS component for BS-RIS channel
d_BS = lambda/2; % BS antenna spacing
N_sqrt = round(sqrt(N)); % RIS as square array
d_RIS = lambda/2; % RIS element spacing

H_LoS = zeros(N, M);
vec_BS_to_RIS = RIS_pos - BS_pos;
distance = norm(vec_BS_to_RIS);

for m = 1:M
    % BS antenna position (ULA along y-axis)
    BS_ant_pos = BS_pos + [0, (m-1)*d_BS, 0];
    
    for n = 1:N
        % RIS element position (UPA)
        row = mod(n-1, N_sqrt);
        col = floor((n-1)/N_sqrt);
        RIS_ant_pos = RIS_pos + [row*d_RIS, col*d_RIS, 0];
        
        % Distance for this antenna pair
        ant_distance = norm(RIS_ant_pos - BS_ant_pos);
        phase = 2*pi*ant_distance/lambda;
        H_LoS(n, m) = exp(-1j*phase);
    end
end

H_LoS = H_LoS / sqrt(N*M); % Normalize
end