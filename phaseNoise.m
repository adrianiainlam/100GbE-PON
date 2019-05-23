function [xPN, phasenoise]  = phaseNoise(x, linewidthTx, linewidthLO, Tsamp)
  %% Simulates laser phase noise.
  %% Inputs:
  %%  - x: input waveform
  %%  - linewidthTx: Tx laser linewidth (Hz)
  %%  - linewidthLO: Rx LO laser linewidth (Hz)
  %%  - Tsamp: Sampling period (s)
  %% Outputs:
  %%  - xPN: output waveform
  %%  - phasenoise: the actual phase noise added (rad)
  dphiTx = sqrt(2 * pi * linewidthTx * Tsamp) * randn(length(x), 1);
  dphiLO = sqrt(2 * pi * linewidthLO * Tsamp) * randn(length(x), 1);
  phiTx = cumsum(dphiTx);
  phiLO = cumsum(dphiLO);

  phasenoise = phiTx - phiLO;
  xPN = x .* exp(-1j * phasenoise);
end
