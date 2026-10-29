function crn = make_crn(K, N, S)
% Common random numbers used by the design-time distortion estimate (model M2).
crn.s = (randn(K,S) + 1j*randn(K,S))/sqrt(2);
crn.v = (randn(N,S) + 1j*randn(N,S))/sqrt(2);
end
