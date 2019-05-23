function [rc, phiests] = phaseNoiseCorr(r, M, phoffset, blocksize)
  %% Phase noise correction.
  %% Inputs:
  %%  - r: Received symbols
  %%  - M: Order of constellation, M=4 for QPSK
  %%  - phoffset: Phase of the 0th symbol in the constellation
  %%  - blocksize: Viterbi-Viterbi algorithm block size
  %% Outputs:
  %%  - rc: Symbols with corrected phase
  %%  - phiests: Estimates of the phase noise

  n = length(r);
  phiests = zeros(1, n);
  rc = zeros(1, n);

  for l = 1 : blocksize : n
    block = r(l : min(l + blocksize - 1, n));

    sum_M = sum(block .^ M);
    phi_est = angle(sum_M .* exp(1j * M * phoffset)) / M;

    if l > 1
      %% Phase unwrapping
      phi_prev = phiests(l - 1);
      m = floor(0.5 + (phi_prev - phi_est) * M / (2 * pi));
      phi_est = phi_est + m * 2 * pi / M;
    end

    block = block .* exp(1j * -phi_est);
    rc(l : min(l + blocksize - 1, n)) = block;
    phiests(l : min(l + blocksize - 1, n)) = phi_est;
  end
end
