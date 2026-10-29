function [sinr, R_sum, P_BS, P_RIS] = compute_metrics2(h, G, f, Theta, P, w, sigma_v_sq, sigma_sq)
% COMPUTE_METRICS2  SINR, sum-rate and power consumption, equations (2)-(4).
%
%   *** OPTIMISED -- numerically identical to the previous version ***
%
%   The previous implementation was O(N^3) per call, dominated by
%
%       T_ris = P * Theta;
%
%   P and Theta are diagonal, but diag() builds them as FULL N x N matrices,
%   so this line performed a dense 256x256 by 256x256 multiply: about 33
%   MFLOPs to compute what is really an element-wise product of two length-N
%   vectors. Three further lines inside the user loop had the same problem:
%
%       f_k' * T_ris * G        dense 1xN by NxN            -> O(N^2)
%       T_ris * G * w(:,k)      dense NxN by NxM, per user  -> O(N^2 M)
%       norm(T_ris,'fro')       over N^2 entries            -> O(N^2)
%
%   Across a full sweep this function is called on the order of 10^7 times,
%   which is why a run took over an hour.
%
%   This version carries the reflection coefficients as a length-N vector
%   t = diag(P).*diag(Theta), which is the diagonal of Psi, and replaces
%   every dense product by the equivalent element-wise operation. G*w is
%   also hoisted out of the user loop, since it does not depend on k.
%
%   Complexity per call drops from O(N^3) to O(NMK).

    if iscell(w),     w     = cell2mat(w);     end
    if iscell(P),     P     = cell2mat(P);     end
    if iscell(Theta), Theta = cell2mat(Theta); end
    if iscell(h),     h     = cell2mat(h);     end
    if iscell(G),     G     = cell2mat(G);     end
    if iscell(f),     f     = cell2mat(f);     end

    K = size(h, 1);

    % Diagonal of Psi = P*Theta, as a vector. Equivalent to diag(P*Theta).
    t = diag(P) .* diag(Theta);              % N x 1

    % Hoisted out of the user loop: independent of k.
    GW  = G * w;                             % N x K
    tGW = bsxfun(@times, t, GW);             % N x K   == T_ris*G*w

    abs_t_sq = abs(t).^2;                    % N x 1

    sinr      = zeros(K, 1);
    RIS_terms = 0;

    for k = 1:K
        % NOTE: both use the conjugate transpose, exactly as in the original,
        % so that h_k' and f_k' below recover the UNCONJUGATED rows h(k,:)
        % and f(k,:). Using .' here instead silently conjugates the channel
        % and changes the SINR.
        h_k = h(k, :)';                      % M x 1
        f_k = f(k, :)';                      % N x 1

        % equiv_channel = h_k' + f_k' * T_ris * G
        equiv_channel = h_k' + (conj(f_k) .* t).' * G;      % 1 x M

        gains  = equiv_channel * w;                          % 1 x K
        pw     = abs(gains).^2;

        signal       = pw(k);
        interference = sum(pw) - signal;

        % norm(f_k' * T_ris)^2 == sum |f_k|^2 |t|^2
        noise = sum(abs(f_k).^2 .* abs_t_sq) * sigma_v_sq + sigma_sq;

        sinr(k) = signal / (interference + noise);

        % norm(T_ris*G*w(:,k))^2
        RIS_terms = RIS_terms + sum(abs(tGW(:, k)).^2);
    end

    % norm(T_ris,'fro')^2 == sum |t|^2
    P_RIS = RIS_terms + sum(abs_t_sq) * sigma_v_sq;
    P_BS  = sum(abs(w(:)).^2);
    R_sum = sum(log2(1 + sinr));
end
