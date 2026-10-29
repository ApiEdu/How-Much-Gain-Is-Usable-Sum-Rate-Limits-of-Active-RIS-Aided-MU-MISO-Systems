function [h, f, G] = Channel_generate_geometric(K, N, M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c, BS_position, RIS_position, user_positions)
% Enhanced channel generation with geometric modeling and corrected path loss

% Constants
c = 3e8;                  % Speed of light
f_c_Hz = f_c * 1e9;       % Frequency in Hz
lambda = c / f_c_Hz;      % Wavelength
d0 = 1;                   % Reference distance
L0_dB = -30;              % Reference path loss in dB
kappa = 1;                % Rician factor

% Convert reference loss to linear
L0_linear = 10^(L0_dB/10);

% Initialize channels
h = zeros(M, K);
f = zeros(N, K);
G = zeros(N, M);

%% ---- BS-RIS channel G ----
path_loss_G_dB = L0_dB - 10*large_fading_DI*log10(Dis_BStoRIS/d0);
path_loss_G_linear = 10^(path_loss_G_dB/10);

H_LoS_G = geometric_LoS_G(M, N, BS_position, RIS_position, lambda);
H_NLoS_G = (randn(N, M) + 1j*randn(N, M)) / sqrt(2);

G = sqrt(path_loss_G_linear) * ...
    (sqrt(kappa/(kappa+1)) * H_LoS_G + sqrt(1/(kappa+1)) * H_NLoS_G);

%% ---- User channels ----
for k = 1:K
    % --- BS-User channel h_k ---
    path_loss_h_dB = L0_dB - 10*large_fading_AI*log10(Dis_BStoUser(k)/d0);
    path_loss_h_linear = 10^(path_loss_h_dB/10);
    
    H_LoS_h = geometric_LoS_h(M, BS_position, user_positions(k,:), lambda);
    H_NLoS_h = (randn(M, 1) + 1j*randn(M, 1)) / sqrt(2);
    
    h(:, k) = sqrt(path_loss_h_linear) * ...
        (sqrt(kappa/(kappa+1)) * H_LoS_h + sqrt(1/(kappa+1)) * H_NLoS_h);
    
    % --- RIS-User channel f_k ---
    path_loss_f_dB = L0_dB - 10*large_fading_DI*log10(Dis_RIStoUser(k)/d0);
    path_loss_f_linear = 10^(path_loss_f_dB/10);
    
    H_LoS_f = geometric_LoS_f(N, RIS_position, user_positions(k,:), lambda);
    H_NLoS_f = (randn(N, 1) + 1j*randn(N, 1)) / sqrt(2);
    
    f(:, k) = sqrt(path_loss_f_linear) * ...
        (sqrt(kappa/(kappa+1)) * H_LoS_f + sqrt(1/(kappa+1)) * H_NLoS_f);
end

end
